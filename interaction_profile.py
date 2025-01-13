import os
import numpy as np
import matplotlib.pyplot as plt

class InteractionProfilePlotter:
    def __init__(self, scores_folder, plots_folder):
        self.scores_folder = scores_folder
        self.plots_folder = plots_folder
        os.makedirs(self.plots_folder, exist_ok=True)

    def read_scores(self, scores_file):
        """Read scores from a file."""
        with open(scores_file, 'r') as f:
            scores = [float(line.strip()) for line in f.readlines()]
        return scores

    def plot_profile(self, base_pair, scores):
        """Plot interaction profile for a base pair."""
        distances = np.linspace(0, 20, len(scores))

        plt.figure(figsize=(8, 6))
        plt.plot(distances, scores, color="blue", linewidth=2, label=f"{base_pair}")
        plt.title(f"Interaction Profile for {base_pair}", fontsize=14)
        plt.xlabel("Distance (Å)", fontsize=12)
        plt.ylabel("Score", fontsize=12)
        plt.grid(alpha=0.4)
        plt.legend()

        output_path = os.path.join(self.plots_folder, f"{base_pair}_profile.png")
        plt.savefig(output_path)
        plt.close()

    def generate_plots(self):
        """Generate plots for all base pairs."""
        for file_name in os.listdir(self.scores_folder):
            if file_name.endswith("_scores.txt"):
                base_pair = file_name.replace("_scores.txt", "")
                scores_file = os.path.join(self.scores_folder, file_name)
                scores = self.read_scores(scores_file)
                self.plot_profile(base_pair, scores)

# Main function to create plots
def main():
    scores_folder = "scores"
    plots_folder = "plots"

    plotter = InteractionProfilePlotter(scores_folder, plots_folder)
    plotter.generate_plots()
    print(f"Plots saved in the '{plots_folder}' folder.")

if __name__ == "__main__":
    main()
