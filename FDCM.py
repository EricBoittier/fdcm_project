import argparse
import configparser
import os
import sys

from analyse_scan import analyse, get_path_neighbours
from job_maker import template_concerted, template_scan, template_morton, template_fit


def template_neighbours(args):
    paths, neighbours = get_path_neighbours(args)

    if not os.path.exists(args.job_folder):
        os.makedirs(args.job_folder)

    n_jobs = len(paths)

    for i, (path, neighbour) in enumerate(zip(paths, neighbours)):
        if i == 0:
            # print(i, path, neighbour)
            tmp_str = template_fit(args, paths[i], paths[i+1], first=True)
            f = open(os.path.join(args.job_folder, f"frame_{path[i]}_{path[i+1]}.sh"))
            f.write(tmp_str)

            for n in neighbour:
                tmp_str = template_fit(args, paths[i], n)
                f_.write()


        elif i < n_jobs:
            pass
        else:
            pass


def main(argv=None):
    # Do argv default this way, as doing it in the functional
    # declaration sets it at compile time.
    if argv is None:
        argv = sys.argv

    # Parse any conf_file specification
    # We make this parser with add_help=False so that
    # it doesn't parse -h and print help.
    conf_parser = argparse.ArgumentParser(
        description=__doc__,  # printed with -h/--help
        # Don't mess with format of description
        formatter_class=argparse.RawDescriptionHelpFormatter,
        # Turn off help, so we print all options in response to -h
        add_help=False
    )
    conf_parser.add_argument("-c", "--conf_file",
                             help="Specify config file", metavar="FILE")
    args, remaining_argv = conf_parser.parse_known_args()

    defaults = {"option": "default"}

    if args.conf_file:
        config = configparser.ConfigParser()
        config.read([args.conf_file])
        defaults.update(dict(config.items("Defaults")))

    # Parse rest of arguments
    # Don't suppress add_help here so it will handle -h
    parser = argparse.ArgumentParser(
        # Inherit options from config_parser
        parents=[conf_parser]
    )
    parser.set_defaults(**defaults)
    parser.add_argument("--option")
    args = parser.parse_args(remaining_argv)

    if args.job_type == "scan":
        output = template_concerted(args)
    elif args.job_type == "concerted":
        output = template_scan(args)
    elif args.job_type == "morton":
        output = template_morton(args)
    elif args.job_type == "neighbours":
        output = template_neighbours(args)
    elif args.job_type == "analysis":
        analyse(args)
    else:
        print("No job type specified")


if __name__ == "__main__":
    main()
