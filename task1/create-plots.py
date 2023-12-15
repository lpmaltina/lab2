from typing import Dict
import matplotlib.pyplot as plt

DIMS = (512, 1024, 2048, 4096, 8192)
MODES = ("rows", "clms", "blocks")
DELIMITER = ","


def read_timings(timings, filename: str) -> Dict[str, Dict[str, Dict[int, float]]]:
    mode = filename.split(".")[0].split("-")[-1]
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                rows, cols, processes, time = line.split(", ")
                rows = int(rows)
                cols = int(cols)
                processes = int(processes)
                time = float(time)
                key = f"{rows}-{cols}-{mode}"
                if key in timings:
                    timings[key]["time"][processes] = time
                else:
                    timings[key] = {"time": {processes: time}}
    return timings


def count_speedup_and_efficiency(
        timings: Dict[str, Dict[str, Dict[int, float]]]
        ):
    for dim in DIMS:
        for mode in MODES:
            key = f"{dim}-{dim}-{mode}"
            one_process_time = timings[key]["time"][1]
            timings[key]["speedup"] = {}
            timings[key]["efficiency"] = {}
            for threads, info in timings[key]["time"].items():
                timings[key]["speedup"][threads] = one_process_time / info
                timings[key]["efficiency"][threads] = (
                    timings[key]["speedup"][threads] / threads
                )


def create_plot(
        timings: Dict[str, Dict[str, Dict[int, float]]],
        feature: str
    ):
    _, axs = plt.subplots(
        1,
        len(MODES),
        figsize=(5 * len(MODES), 5)
    )
    plt.suptitle(feature.title())
    for dim in DIMS:
        for i, mode in enumerate(MODES):
            key = f"{dim}-{dim}-{mode}"
            label = f"{dim}x{dim}"
            axs[i].set_title(mode)
            axs[i].plot(
                timings[key][feature].keys(),
                timings[key][feature].values(),
                label=label
            )
    plt.legend()
    plt.savefig(f"plots/{feature}.png")


timings = {}
for mode in MODES:
    read_timings(timings, f"timings/timings-{mode}.txt")
count_speedup_and_efficiency(timings)

for feature in ("speedup", "efficiency"):
    create_plot(timings, feature)