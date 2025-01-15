import os
import numpy as np

class RNAGibbsFreeEnergyCalculator:
    def __init__(self, scores_folder, distance_threshold=20, sequence_threshold=4):
        self.scores_folder = scores_folder
        self.distance_threshold = distance_threshold
        self.sequence_threshold = sequence_threshold
        self.distance_bins = np.linspace(0, distance_threshold, 21)
        self.scores = self.load_scores()

    def load_scores(self):
        """Load precomputed scores for all base pairs."""
        scores = {}
        for file_name in os.listdir(self.scores_folder):
            if file_name.endswith("_scores.txt"):
                base_pair = file_name.replace("_scores.txt", "")
                with open(os.path.join(self.scores_folder, file_name), 'r') as f:
                    scores[base_pair] = np.array([float(line.strip()) for line in f.readlines()])
        return scores

    def interpolate_score(self, base_pair, distance):
        """Perform linear interpolation for the score."""
        if base_pair not in self.scores:
            return 10  # Arbitrary high score for missing base pairs

        bin_edges = self.distance_bins
        bin_scores = self.scores[base_pair]
        
        if distance < bin_edges[0] or distance > bin_edges[-1]:
            return 10  # Maximum score for out-of-range distances

        bin_index = np.digitize(distance, bin_edges) - 1
        if bin_index >= len(bin_scores) - 1:
            return bin_scores[-1]

        # Linear interpolation
        r1, r2 = bin_edges[bin_index], bin_edges[bin_index + 1]
        s1, s2 = bin_scores[bin_index], bin_scores[bin_index + 1]
        return s1 + (distance - r1) * (s2 - s1) / (r2 - r1)

    def parse_pdb(self, pdb_file):
        """Parse a PDB file and extract C3' atom coordinates."""
        atoms = []
        with open(pdb_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM") and " C3' " in line:
                    res_name = line[17:20].strip()
                    coord = (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54])
                    )
                    atoms.append((res_name, coord))
        return atoms

    def compute_gibbs_free_energy(self, pdb_file):
        """Compute the normalized Gibbs free energy for the RNA structure."""
        atoms = self.parse_pdb(pdb_file)
        total_score = 0
        valid_pairs = 0

        for i, (res1, coord1) in enumerate(atoms):
            for j, (res2, coord2) in enumerate(atoms):
                if abs(i - j) < self.sequence_threshold:
                    continue

                distance = np.linalg.norm(np.array(coord1) - np.array(coord2))
                if distance > self.distance_threshold:
                    continue

                base_pair = res1[0] + res2[0]
                score = self.interpolate_score(base_pair, distance)
                total_score += score
                valid_pairs += 1

        return total_score / valid_pairs if valid_pairs > 0 else float('inf')

# Main function to compute Gibbs free energy for multiple structures
def main():
    scores_folder = "scores"
    structures_folder = "rna_structures"  # Folder containing the predicted RNA structures
    output_file = "gibbs_free_energy_results.txt"

    calculator = RNAGibbsFreeEnergyCalculator(scores_folder)
    results = []

    for pdb_file in os.listdir(structures_folder):
        if pdb_file.endswith(".pdb"):
            pdb_path = os.path.join(structures_folder, pdb_file)
            gibbs_free_energy = calculator.compute_gibbs_free_energy(pdb_path)
            results.append((pdb_file, gibbs_free_energy))

    # Sort results by Gibbs free energy (ascending order, lower is better)
    results.sort(key=lambda x: x[1])

    # Write results to a file
    with open(output_file, 'w') as f:
        f.write("Normalized Gibbs Free Energy Results:\n")
        for pdb_file, energy in results:
            f.write(f"{pdb_file}: {energy}\n")

    print(f"Results saved to '{output_file}'.")
    print("Normalized Gibbs Free Energy of all structures:")
    for pdb_file, energy in results:
        print(f"{pdb_file}: {energy}")

    print("\nTop structure:")
    print(f"{results[0][0]} with Gibbs free energy: {results[0][1]}")

if __name__ == "__main__":
    main()
