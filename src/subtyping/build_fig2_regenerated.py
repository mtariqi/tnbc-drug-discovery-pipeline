import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

NAVY = "#1F3864"
GOLD = "#8A6D00"
TEAL = "#2E7D6E"
GREY = "#8C8C8C"

fig, ax = plt.subplots(figsize=(9, 5.5))

stages = [
    "Sequenced\nTCGA-BRCA cases",
    "PAM50 basal-like\n(BRCA_Subtype_PAM50)",
    "+ real mutation\ndata available",
    "\u2265 2 candidate\ndrug-mapped genes",
]
values = [1087, 192, 168, 18]
colors = [GREY, GOLD, TEAL, NAVY]

y_pos = np.arange(len(stages))[::-1]
bars = ax.barh(y_pos, values, color=colors, height=0.6)

for y, v in zip(y_pos, values):
    ax.text(v + max(values) * 0.015, y, f"{v}", va="center", ha="left", fontsize=13, fontweight="bold")

ax.set_yticks(y_pos)
ax.set_yticklabels(stages, fontsize=11)
ax.set_xlabel("Number of participants", fontsize=11)
ax.set_title("Cohort-scale screening funnel for TNBC combination-therapy candidates", fontsize=13, pad=15)
ax.set_xlim(0, max(values) * 1.18)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("fig2_regenerated.png", dpi=200, bbox_inches="tight")
print("Saved.")

# Overlap audit
import matplotlib.transforms as mtransforms
print("Real values plotted:", values)
