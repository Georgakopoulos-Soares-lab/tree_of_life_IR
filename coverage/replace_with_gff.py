import json
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--schedule", type=str)
    args = parser.parse_args()
    schedule = args.schedule
    with open(schedule, "r") as f:
        data = json.load(f)
    new_schedule = defaultdict(list)
    for bucket_id, bucket in tqdm(data.items()):
        for file in bucket:
            file_gff = file.replace("fna", "gff")
            if Path(file_gff).is_file():
                new_schedule[bucket_id].append(file_gff)
    new_schedule = dict(new_schedule)
    with open(schedule.replace(".json", ".rescheduled.json"), "w") as f:
        json.dump(new_schedule, f, indent=4)



    
