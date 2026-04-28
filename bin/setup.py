#!/user/bin/python3



from utils import *
import argparse



# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--organism', type=str, default='homo_sapiens', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    return parser.parse_args()

# Parse command line arguments
args = parse_arguments()

# Read the JSON tree file
#with open(args.tree_file, 'r') as file:
   # tree = json.load(file)

# Set organism and census_version from arguments
organism = args.organism
census_version = args.census_version

# Get model file link and download
model_path = setup(organism=organism, version=census_version)
