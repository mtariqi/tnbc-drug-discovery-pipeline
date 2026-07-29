"""
NOT PART OF THE PACKAGE. Verification-only harness.

pytest is not installed in this sandbox and there is no network to
install it. This script implements just enough of pytest's fixture
mechanism (fixture registration/resolution, monkeypatch, tmp_path,
raises) to actually execute the real tests/test_*.py files unmodified,
so their assertions are genuinely checked rather than assumed correct
from reading the code. Run `pytest tests/` for real once pytest is
available in your environment -- this harness is a stand-in, not a
replacement.
"""

import importlib.util
import inspect
import shutil
import sys
import tempfile
import types
from contextlib import contextmanager
from pathlib import Path


# ---------------------------------------------------------------------
# Minimal fake `pytest` module
# ---------------------------------------------------------------------

fake_pytest = types.ModuleType("pytest")


def fixture(func=None, **kwargs):
    def decorator(f):
        f._is_fixture = True
        return f
    if func is not None:
        return decorator(func)
    return decorator


@contextmanager
def raises(exc_type, match=None):
    import re
    try:
        yield
    except exc_type as e:
        if match and not re.search(match, str(e)):
            raise AssertionError(f"exception message {str(e)!r} did not match pattern {match!r}")
        return
    else:
        raise AssertionError(f"expected {exc_type.__name__} to be raised, nothing was")


fake_pytest.fixture = fixture
fake_pytest.raises = raises
sys.modules["pytest"] = fake_pytest


# ---------------------------------------------------------------------
# Minimal monkeypatch + tmp_path builtin-fixture stand-ins
# ---------------------------------------------------------------------

class FakeMonkeyPatch:
    def __init__(self):
        self._undo = []
        self._env_undo = []

    def setattr(self, target, name, value=None):
        if value is None:
            # target is "obj.attr" as a single dotted string form not used here
            raise NotImplementedError
        original = getattr(target, name)
        self._undo.append((target, name, original))
        setattr(target, name, value)

    def setenv(self, name, value):
        import os
        had_original = name in os.environ
        original = os.environ.get(name)
        self._env_undo.append((name, had_original, original))
        os.environ[name] = value

    def undo(self):
        for target, name, original in reversed(self._undo):
            setattr(target, name, original)
        self._undo.clear()
        import os
        for name, had_original, original in reversed(self._env_undo):
            if had_original:
                os.environ[name] = original
            else:
                os.environ.pop(name, None)
        self._env_undo.clear()


_tmp_dirs = []


def make_tmp_path():
    d = tempfile.mkdtemp()
    _tmp_dirs.append(d)
    return Path(d)


# ---------------------------------------------------------------------
# Fixture resolution + test execution
# ---------------------------------------------------------------------

def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def resolve_and_call(func, fixture_sources, cache):
    """
    fixture_sources: a list of objects searched IN ORDER for a fixture matching a
    requested parameter name -- e.g. [test_class_instance, test_module, conftest_module].
    Real pytest resolves fixtures with this kind of scoping (test-local before conftest);
    this harness's original version only ever checked a single conftest module, which
    silently can't find a fixture defined inside the test file itself (the actual,
    real situation in test_redundancy_analyzer.py's module-scoped `analyzer` fixture).
    """
    sig = inspect.signature(func)
    kwargs = {}
    monkeypatch_instance = None
    for pname in sig.parameters:
        if pname == "self":
            continue
        if pname == "monkeypatch":
            monkeypatch_instance = FakeMonkeyPatch()
            kwargs["monkeypatch"] = monkeypatch_instance
        elif pname == "tmp_path":
            kwargs["tmp_path"] = make_tmp_path()
        elif pname in cache:
            kwargs[pname] = cache[pname]
        else:
            fixture_func = None
            for source in fixture_sources:
                if source is not None and hasattr(source, pname):
                    candidate = getattr(source, pname)
                    if callable(candidate) and getattr(candidate, "_is_fixture", False):
                        fixture_func = candidate
                        break
            if fixture_func is None:
                raise RuntimeError(f"no fixture found for parameter {pname!r} needed by {func.__name__}")
            value = resolve_and_call(fixture_func, fixture_sources, cache)
            cache[pname] = value
            kwargs[pname] = value
    try:
        return func(**kwargs)
    finally:
        if monkeypatch_instance is not None:
            monkeypatch_instance.undo()


def run_test_file(test_path, conftest):
    """
    Discovers both flat test_* functions AND class-based tests
    (class Test...: def test_...(self, fixture_name)), since the real test suite here
    uses the class-based style with a fixture defined inside the same module -- the
    original harness only handled the flat-function case and would have silently found
    zero tests in test_redundancy_analyzer.py, not actually verified anything in it.
    """
    module = load_module(test_path, test_path.stem)
    results = []

    # Flat, module-level test functions
    for name, obj in vars(module).items():
        if name.startswith("test_") and callable(obj) and inspect.isfunction(obj):
            cache = {}
            try:
                resolve_and_call(obj, [module, conftest], cache)
                results.append((name, "PASS", None))
            except Exception as e:
                results.append((name, "FAIL", e))

    # Class-based tests (unittest/pytest-style `class Test...`)
    for cls_name, cls in vars(module).items():
        if inspect.isclass(cls) and cls_name.lower().startswith("test"):
            instance = cls()
            for method_name in dir(instance):
                if method_name.startswith("test_"):
                    method = getattr(instance, method_name)
                    cache = {}
                    try:
                        resolve_and_call(method, [instance, module, conftest], cache)
                        results.append((f"{cls_name}.{method_name}", "PASS", None))
                    except Exception as e:
                        results.append((f"{cls_name}.{method_name}", "FAIL", e))

    return results


def main():
    repo_root = Path(__file__).parent
    sys.path.insert(0, str(repo_root))

    conftest_path = repo_root / "tests" / "conftest.py"
    conftest = load_module(conftest_path, "conftest") if conftest_path.exists() else None
    if conftest is None:
        print("Note: no tests/conftest.py found -- proceeding without one (fixtures defined "
              "inside individual test modules are still resolved; this is only a problem if "
              "a test relies on a *shared*, cross-file fixture that genuinely lives in conftest.py).")

    test_files = sorted((repo_root / "tests").glob("test_*.py"))
    total_pass, total_fail = 0, 0
    for test_path in test_files:
        print(f"\n=== {test_path.name} ===")
        results = run_test_file(test_path, conftest)
        for name, status, err in results:
            if status == "PASS":
                print(f"  PASS  {name}")
                total_pass += 1
            else:
                print(f"  FAIL  {name}: {err!r}")
                import traceback
                traceback.print_exception(type(err), err, err.__traceback__)
                total_fail += 1

    print(f"\n{total_pass} passed, {total_fail} failed")

    for d in _tmp_dirs:
        shutil.rmtree(d, ignore_errors=True)

    if total_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
