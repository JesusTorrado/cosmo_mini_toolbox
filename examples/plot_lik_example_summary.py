import os
import sys
import matplotlib.pyplot as plt

sys.path.append("../src")
from Chain import Chain
from plot_lik import plot_lik_2D

base_folder = "./chains"
chain = "planck_WP"
chains = [Chain(os.path.join(base_folder, chain))]
params=["H0",     "omega_b"]
labels=[r"$H_0$", r"$\omega_b$"]

fig, axarr = plt.subplots(2,2)
axes_locations = {                    "profile":  axarr[0,1],
                  "mean": axarr[1,0], "marginal": axarr[1,1]}
for mode, axes in axes_locations.items():
    ax, options = plot_lik_2D(mode, chains, params=params, labels=labels, format = "-loglik", dpi=200, fontsize_labels=14, fontsize_ticks=8,
                save=0, axes=axes
                )
    axes.set_title(mode.title(), fontdict={'fontsize':16})

# Text:
text = ("Chain:\n %s\n\n"%[c.name() for c in chains] +
        "Parameters:\n%s"%(params))
axarr[0,0].set_axis_off()
axarr[0,0].text(0, 1, text, weight="bold", verticalalignment="top")

# Plot
plt.tight_layout()
plt.savefig("summary.png", **options)
plt.show()
plt.close()


