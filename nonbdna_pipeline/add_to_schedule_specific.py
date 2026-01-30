def main():
    import argparse
    import json
    from pathlib import Path
    from termcolor import colored
    from scheduling import MiniBucketScheduler
    parser = argparse.ArgumentParser(description="Rescheduler")
    parser.add_argument("schedule", type=str)
    parser.add_argument("-i", "--indir", type=str, default="primates_ref/fasta/")
    targets = [149, 192, 272, 59]
    args = parser.parse_args()
    schedule = Path(args.schedule).resolve()
    if not schedule.is_file():
        raise FileNotFoundError(f"Schedule `{schedule}` doesn't exist.")
    indir = Path(args.indir).resolve()
    if not indir.is_dir():
        raise FileNotFoundError(f"Directory `{indir}` doesn't exist.")
    scheduler = MiniBucketScheduler()
    with open(schedule, mode="r", encoding="UTF-8") as f:
        schedule_data = json.load(f)

    infiles = [infile for infile in indir.glob("*.fna")]
    infiles += infiles + [infile for infile in indir.glob("*.fa")]
    infiles += infiles + [infile for infile in indir.glob("*.fna.gz")]
    infiles += infiles + [infile for infile in indir.glob("*.fa.gz")]
    infiles = list(map(str, infiles))
    total_infiles = len(infiles)
    assign_to_bucket = dict()
    # Assign new files to random buckets proportionally
    targets = list(map(str, targets))
    for i, infile in enumerate(infiles):
        assign_to_bucket.setdefault(targets[i%len(targets)], []).append(infile)
    print(assign_to_bucket)

    print(colored(f"Total infiles detected: {total_infiles}.", "yellow"))
    # Append new files
    for bucket_id, bucket in assign_to_bucket.items():
        for file in bucket:
            schedule_data[bucket_id].append(file)

    dest = schedule.parent.joinpath(schedule.name.replace(".json", f".updated_{total_infiles}.json"))
    with open(dest, mode="w", encoding="UTF-8") as fout:
        json.dump(schedule_data, fout, indent=4)

    print(colored(f"Updated schedule at {dest}!", "green"))
if __name__ == "__main__": main()
