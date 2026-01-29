import json 
import uuid
from pathlib import Path 
from termcolor import colored
from scheduling import MiniBucketScheduler
def main():
    import argparse; parser = argparse.ArgumentParser()
    parser.add_argument("schedule", type=str)
    parser.add_argument("--total_buckets", type=int, default=10)
    args = parser.parse_args()
    schedule = Path(args.schedule).resolve()
    if not schedule.is_file():
        raise ValueError(f"Missing schedule file `{schedule}`.")
    
    with open(schedule, mode="r", encoding="UTF-8") as f:
        schedule_data = json.load(f)
    scheduled_files = dict()
    extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
    for _, bucket in schedule_data.items():
        for file in bucket:
            scheduled_files[extract_id(file)] = file
    infiles = [infile for infile in Path().cwd().glob("design_bucket_*.csv") if infile.is_file()]
    print(colored(f"Total infiles detected: {len(infiles)}.", "blue"))
    for infile in infiles:
        with open(infile, mode="r", encoding="UTF-8") as fin:
            for line in fin:
                processed_accession_id = extract_id(line.strip().split(",")[0])
                scheduled_files.pop(processed_accession_id)
    # rescheduled_file = schedule.parent.joinpath(schedule.name.replace(".json", "_rescheduled.json"))
    print(colored(f"Remaining files: {len(scheduled_files)} to be rescheduled.", "blue"))
    scheduler = MiniBucketScheduler()
    scheduled = scheduler.schedule(files=list(scheduled_files.values()), total_buckets=args.total_buckets)
    unique_id = str(uuid.uuid4()).replace("-", "")
    destination = schedule.parent.joinpath(f"rescheduled_schedule_{args.total_buckets}_{unique_id}.json")
    scheduler.saveas(scheduled, destination)
    print(colored(f"New schedule has been saved at {destination}", "green"))
    print(colored("DONE!", "green"))

if __name__ == "__main__": main()