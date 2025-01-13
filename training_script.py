import os
import numpy as np

# Constants
DISTANCE_THRESHOLD = 20  # Maximum distance in Å
SEQUENCE_THRESHOLD = 4   # Minimum sequence separation
BASE_PAIRS = ["AA", "AU", "AC", "AG", "UU", "UC", "UG", "CC", "CG", "GG"]
DISTANCE_BINS = np.linspace(0, DISTANCE_THRESHOLD, 21)  # 20 intervals

class RNAInteractionScorer:
    def __init__(self, distance_threshold, sequence_threshold, base_pairs, distance_bins):
        self.distance_threshold = distance_threshold
        self.sequence_threshold = sequence_threshold
        self.base_pairs = base_pairs
        self.distance_bins = distance_bins
        self.num_bins = len(distance_bins) - 1

    @staticmethod
    def compute_distance(coord1, coord2):
        """Compute Euclidean distance between two coordinates."""
        return np.linalg.norm(np.array(coord1) - np.array(coord2))

    def parse_pdb(self, file_path):
        """Parse a PDB file and extract C3' atom coordinates."""
        atoms = []
        with open(file_path, 'r') as file:
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

    def count_distances(self, atoms):
        """Count distances for base pairs and reference pairs."""
        counts = {bp: np.zeros(self.num_bins) for bp in self.base_pairs}
        ref_counts = np.zeros(self.num_bins)

        for i, (res1, coord1) in enumerate(atoms):
            for j, (res2, coord2) in enumerate(atoms):
                if abs(i - j) < self.sequence_threshold:
                    continue

                distance = self.compute_distance(coord1, coord2)
                if distance > self.distance_threshold:
                    continue

                distance_bin = np.digitize(distance, self.distance_bins) - 1
                if distance_bin >= self.num_bins:
                    continue

                base_pair = res1[0] + res2[0]
                ref_counts[distance_bin] += 1

                if base_pair in counts:
                    counts[base_pair][distance_bin] += 1

        return counts, ref_counts

    def compute_scores(self, counts, ref_counts):
        """Compute log-ratio scores."""
        scores = {}
        total_ref = ref_counts.sum()

        for bp, bp_counts in counts.items():
            scores[bp] = []
            total_bp = bp_counts.sum()

            for i in range(self.num_bins):
                f_obs = bp_counts[i] / total_bp if total_bp > 0 else 0
                f_ref = ref_counts[i] / total_ref if total_ref > 0 else 0
                score = -np.log(f_obs / f_ref) if f_obs > 0 and f_ref > 0 else 10
                scores[bp].append(min(score, 10))

        return scores

    def save_scores(self, scores, output_dir):
        """Save scoring profiles to files."""
        os.makedirs(output_dir, exist_ok=True)
        for bp, bp_scores in scores.items():
            file_path = os.path.join(output_dir, f"{bp}_scores.txt")
            with open(file_path, 'w') as f:
                for score in bp_scores:
                    f.write(f"{score}\n")

# Main function
def main():
    pdb_folder = "pdb_files"
    output_folder = "scores"

    scorer = RNAInteractionScorer(
        distance_threshold=DISTANCE_THRESHOLD,
        sequence_threshold=SEQUENCE_THRESHOLD,
        base_pairs=BASE_PAIRS,
        distance_bins=DISTANCE_BINS
    )

    all_counts = {bp: np.zeros(len(DISTANCE_BINS) - 1) for bp in BASE_PAIRS}
    all_ref_counts = np.zeros(len(DISTANCE_BINS) - 1)

    for pdb_file in os.listdir(pdb_folder):
        if pdb_file.endswith(".pdb"):
            file_path = os.path.join(pdb_folder, pdb_file)
            atoms = scorer.parse_pdb(file_path)
            counts, ref_counts = scorer.count_distances(atoms)

            for bp in BASE_PAIRS:
                all_counts[bp] += counts[bp]
            all_ref_counts += ref_counts

    scores = scorer.compute_scores(all_counts, all_ref_counts)
    scorer.save_scores(scores, output_folder)

if __name__ == "__main__":
    main()
