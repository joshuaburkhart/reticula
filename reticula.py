import argparse
import os
import sys
from pathlib import Path
from string import Template

from dotenv import load_dotenv


def execute_program() -> None:
    os.chdir(os.getenv(varname='buildPath'))
    os.system(f"sh {os.getenv(varname='buildInstFilNme')}")


def generate_program_file_from_template(template_path: str,
                                        instance_name: str,
                                        instance_permissions: int) -> None:
    with open(template_path, 'r') as build_template_filepointer:
        build_template = build_template_filepointer.read()
        build_script = Template(build_template).substitute(os.environ)
    instance_path = f"{os.getenv(varname='buildPath')}/{instance_name}"
    with open(instance_path, 'w') as build_script_fp:
        build_script_fp.write(build_script)
        os.chmod(instance_path, instance_permissions)


def generate_program_files() -> None:
    generate_program_file_from_template(os.getenv(varname='buildTemplPath'),
                                        os.getenv(varname='buildInstFilNme'),
                                        0o555)  # set permissions to "r-xr-xr-x"
    generate_program_file_from_template(os.getenv(varname='dockerTemplPath'),
                                        os.getenv(varname='dockerInstFilNme'),
                                        0o444)  # set permissions to "r--r--r--"
    generate_program_file_from_template(os.getenv(varname='slurmTemplPath'),
                                        os.getenv(varname='slurmInstFilNme'),
                                        0o555)  # set permissions to "r-xr-xr-x"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configuration", "-c", type=str, default="reticula_config.env",
                        help="Supply the full path to an environment configuration (.env) file.")
    return parser


def main() -> None:
    args = build_parser().parse_args(sys.argv[1:])
    env_path = Path(args.configuration)
    load_dotenv(dotenv_path=env_path)
    generate_program_files()
    execute_program()
