import argparse
import configparser
import sys

from analyse_scan import analyse
from job_maker import template_concerted, template_scan, template_morton, \
    template_neighbours_from_2drmsd, template_neighbours_from_anytree_and_G


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

    #  Banch based on job type
    if args.job_type == "scan":
        output = template_concerted(args)
        print(output)
    elif args.job_type == "concerted":
        output = template_scan(args)
    elif args.job_type == "morton":
        output = template_morton(args)
        print(output)
    elif args.job_type == "graph":
        output = template_neighbours_from_anytree_and_G(args)
    elif args.job_type == "analyse":
        analyse(args)
    else:
        print("No job type specified")


if __name__ == "__main__":
    main()
