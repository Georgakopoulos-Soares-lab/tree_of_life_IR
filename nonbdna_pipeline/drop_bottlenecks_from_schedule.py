def main():
    import json 
    import os
    from termcolor import colored
    from pathlib import Path
    import argparse 
    parser = argparse.ArgumentParser()
    parser.add_argument("schedule", type=str)
    parser.add_argument("--indir", "-i", type=str)
    args = parser.parse_args()
    indir = Path(args.indir).resolve()
    with open(args.schedule, mode="r", encoding="UTF-8") as f:
        schedule_data = json.load(f)
    log_indir = indir.joinpath("log_debug_nonbdna")
    job_status_indir = indir.joinpath("status")
    infiles = [infile for infile in log_indir.glob("mindi_*.log")]
    # Fetch bottleneck
    bottlenecks = set()
    extract_id = lambda x: "_".join(Path(x).name.split("_")[:2])
    for infile in infiles:
        current_bid = infile.name.split("_")[2].split(".")[0]
        bucket_status = job_status_indir.joinpath(f"bucket_{current_bid}_status.completed")
        with open(infile, mode="r", encoding="UTF-8") as f:
            for line in f:
                pass 
            if "Processing accession" in line:
                bottlenecks.add(extract_id(line.split("Processing accession: ")[1].strip().split(" ")[0]))
                if bucket_status.is_file():
                    print(colored(f"Job Alert! Removing falsely reported bucket status file `{bucket_status}`. Reason: Bucket {current_bid} is not complete.", "red"))
                    os.remove(bucket_status)
    total_bottlenecks = len(bottlenecks)
    print(colored(f"Total bottlenecks found: {total_bottlenecks}.", "yellow"))
    if total_bottlenecks == 0:
        print("No bottlenecks detected! Happy?")
        print("Not Happy? Exiting...")
        return

    new_schedule = dict()
    for bucket_id, bucket in schedule_data.items():
        for infile in bucket:
            accession_id = extract_id(infile)
            if accession_id in bottlenecks:
                print(colored(f"Omitting bottleneck {infile}...", "yellow"))
                continue
            new_schedule.setdefault(bucket_id, []).append(infile)
    destination = Path(args.schedule).parent.joinpath(f"dropped_bottlenecks_{Path(args.schedule).name}")
    with open(destination, mode="w", encoding="UTF-8") as f:
        json.dump(new_schedule, f, indent=4)
    print(colored(f"New schedule without bottlenecks saved at {destination}", "green"))
if __name__ == "__main__": main()
