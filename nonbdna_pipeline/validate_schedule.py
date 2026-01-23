import json 
from termcolor import colored 

def validate_schedule(schedule, accession_list):
    filtered_accessions = dict()
    with open(accession_list, mode="r", encoding="UTF-8") as f:
        for line in f:
            filtered_accessions.update({line.strip(): False})
    
    with open(schedule, mode="r", encoding="UTF-8") as f:
        data = json.load(f)

    for bucket, files in data.items():
        for file in files:
            filtered_accessions[file] = True
            
    for accession, found in filtered_accessions.items():
        if not found:
            raise ValueError(f"Accession {accession} from filtered list not found in schedule.")
    
    print(colored(f"All checks PASSED!", "green"))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Validate that all accessions in the filtered list are present in the schedule.")
    parser.add_argument("--schedule", type=str, required=True, help="Path to the schedule JSON file.")
    parser.add_argument("--accession_list", type=str, required=True, help="Path to the filtered accession list file.")

    args = parser.parse_args()
    validate_schedule(args.schedule, args.accession_list)