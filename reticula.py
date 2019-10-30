import sys
import argparse

def execute_program():
    pass

def generate_program_files():
    pass

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--configuration", "-c", required=True,
                        help="Please supply the full path to a reticula_config.json5 configuration file.")
    return parser

def main():
    args = build_parser().parse_args(sys.argv[1:])
    generate_program_files(args)
    execute_program(args)
