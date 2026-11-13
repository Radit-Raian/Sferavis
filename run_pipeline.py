"""Run the whole pipeline.

    python run_pipeline.py                                    # asks for RA, Dec and search radius (deg)
    python run_pipeline.py --ra 180 --dec 0 --radius 0.1      # no questions
    python run_pipeline.py --no-prompt                        # reuse the values already in config.json

Works from any folder (all paths are relative to this repo), e.g.  python /path/to/run_pipeline.py
"""
import argparse
import importlib

import helpers

STEPS = ("step2_download_obs", "step3_detect_sources", "step4_survey_match",
         "step5_list_unlabeled_obj", "step6_link_tracklets", "step7_check_results")


def main():
    ap = argparse.ArgumentParser(description="SPHEREx moving-object pipeline")
    ap.add_argument("--ra"), ap.add_argument("--dec")
    ap.add_argument("--radius", type=float, help="search radius in degrees")
    ap.add_argument("--no-prompt", action="store_true", help="do not ask for RA/Dec/radius")
    args = ap.parse_args()

    helpers.check_dependencies()  # clear message instead of an ImportError traceback

    print("\n=== step1_config ===")
    step1 = importlib.import_module("step1_config")
    step1.main(args.ra, args.dec, args.radius, interactive=False if args.no_prompt else None)
    for name in STEPS:
        print(f"\n=== {name} ===")
        importlib.import_module(name).main()


if __name__ == "__main__":
    main()