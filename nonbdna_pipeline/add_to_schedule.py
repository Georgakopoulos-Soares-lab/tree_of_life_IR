def main():
    import argparse
    import json
    from pathlib import Path
    from termcolor import colored
    from scheduling import MiniBucketScheduler
    parser = argparse.ArgumentParser(description="Rescheduler")
    parser.add_argument("schedule", type=str)
    parser.add_argument("-i", "--indir", type=str, default="primates_ref/fasta/")
    parser.add_argument("--total_buckets", "-t", type=int, default=288)
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

    previous_files = [] 
    for _, bucket in schedule_data.items():
        previous_files += [Path(infile) for infile in bucket]
    infiles = [infile for infile in indir.glob("*.fna")]
    infiles += infiles + [infile for infile in indir.glob("*.fa")]
    infiles += infiles + [infile for infile in indir.glob("*.fna.gz")]
    infiles += infiles + [infile for infile in indir.glob("*.fa.gz")]
    total_infiles = len(infiles)
    print(colored(f"Total infiles detected: {total_infiles}.", "yellow"))
    infiles += previous_files
    print(colored(f"Total infiles combined: {len(infiles)}.", "yellow"))
    print(f"Adding these files randomly...")
    infiles = [str(infile) for infile in infiles]
    new_schedule = scheduler.schedule(infiles, total_buckets=args.total_buckets)
    dest = schedule.parent.joinpath(schedule.name.replace(".json", ".filled.json"))
    scheduler.saveas(new_schedule, dest=dest)
    print(colored("Jobs done!", "green"))
if __name__ == "__main__": main()
