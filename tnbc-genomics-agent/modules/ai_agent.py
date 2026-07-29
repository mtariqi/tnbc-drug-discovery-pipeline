"""
ai_agent.py
AI reasoning layer: calls Claude (claude-sonnet-4-6) to synthesize variant
findings, redundancy data, and drug options into a clinician-readable report.
"""

import json
from typing import Dict


class TNBCReasoningAgent:
    """
    Uses the Claude API to reason over structured bioinformatics data
    and generate a clinician-facing narrative summary.

    This module is designed to be called from the dashboard (JavaScript)
    via the Anthropic API, but also exposes a Python-callable interface
    for use in notebooks or CLI pipelines.
    """

    SYSTEM_PROMPT = """You are an expert oncology bioinformatician specializing in
Triple-Negative Breast Cancer (TNBC) genomics. Your role is to analyze structured
genomic data — RTK/nRTK variant profiles, pathway redundancy scores, and drug
interaction data — and produce a clear, evidence-grounded clinical summary.

Guidelines:
- Lead with the most actionable findings.
- Explicitly address pathway redundancy and bypass risk.
- Suggest rational drug combinations where redundancy is detected.
- Distinguish FDA-approved therapies from investigational options.
- Use precise oncology terminology but remain accessible to a clinical audience.
- Conclude with a concise "Priority Action" section.
- Do NOT fabricate clinical data; reason only from the provided structured input.
"""

    @staticmethod
    def build_user_prompt(vcf_summary: Dict, redundancy_report: Dict, patient_id: str = "PATIENT_1") -> str:
        """
        Constructs the structured prompt sent to Claude.
        Converts Python dicts to clean JSON for reliable parsing.
        """
        vcf_json = json.dumps(vcf_summary, indent=2)
        red_json = json.dumps(redundancy_report, indent=2)

        return f"""## TNBC Genomic Analysis Request — Patient: {patient_id}

### VCF Variant Summary
```json
{vcf_json}
```

### RTK/nRTK Redundancy Analysis
```json
{red_json}
```

Please provide:
1. **Variant Significance Summary** – which RTK/nRTK alterations are most clinically significant and why.
2. **Pathway Redundancy Assessment** – which pathways are most redundantly activated and what this means therapeutically.
3. **Top Drug Targets** – rank the top 3 targetable kinases with supporting rationale.
4. **Combination Therapy Recommendations** – specific drug pairs to address redundancy.
5. **Resistance Risk Factors** – known mechanisms that may limit therapeutic efficacy.
6. **Priority Action** – a 2-3 sentence executive summary for clinical decision-making.
"""

    @staticmethod
    def get_js_fetch_snippet(patient_id: str, vcf_summary: Dict, redundancy_report: Dict) -> str:
        """
        Returns the JavaScript fetch() call used in the dashboard artifact.
        Kept here to keep the AI call logic co-located with the Python agent.
        """
        prompt = TNBCReasoningAgent.build_user_prompt(vcf_summary, redundancy_report, patient_id)
        # Escape for embedding in JS template literal
        prompt_escaped = prompt.replace("`", "\\`").replace("${", "\\${")
        return f"""
const response = await fetch("https://api.anthropic.com/v1/messages", {{
  method: "POST",
  headers: {{ "Content-Type": "application/json" }},
  body: JSON.stringify({{
    model: "claude-sonnet-4-6",
    max_tokens: 1000,
    system: `{TNBCReasoningAgent.SYSTEM_PROMPT.replace(chr(10), "\\n")}`,
    messages: [{{ role: "user", content: `{prompt_escaped}` }}]
  }})
}});
const data = await response.json();
const report = data.content.map(b => b.text || "").join("\\n");
""".strip()
