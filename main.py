import argparse
import os
import sys
from pathlib import Path
from string import Template

from dotenv import load_dotenv


def execute_program() -> None:
    os.chdir(os.getenv(key='buildPath'))
    os.system(f"sh {os.getenv(key='buildInstFilNme')}{os.getenv(key='buildInstExt')}")


def generate_program_file_from_template(template_path: str,
                                        instance_name: str,
                                        instance_extension: str,
                                        instance_permissions: int) -> None:
    with open(template_path, 'r') as build_template_filepointer:
        build_template = Template(build_template_filepointer.read())
        build_script = build_template.safe_substitute(os.environ)
    instance_path = f"{os.getenv(key='buildPath')}/{instance_name}{instance_extension}"
    with open(instance_path, 'w') as build_script_fp:
        build_script_fp.write(build_script)
        os.chmod(instance_path, instance_permissions)


def generate_program_files() -> None:
    generate_program_file_from_template(os.getenv(key='buildTemplPath'),
                                        os.getenv(key='buildInstFilNme'),
                                        os.getenv(key='buildInstExt'),
                                        0o744)  # set permissions to "rwxr--r--"
    generate_program_file_from_template(os.getenv(key='dockerTemplPath'),
                                        os.getenv(key='dockerInstFilNme'),
                                        os.getenv(key='dockerInstExt'),
                                        0o644)  # set permissions to "rw-r--r--"
    generate_program_file_from_template(os.getenv(key='slurmTemplPath'),
                                        os.getenv(key='slurmInstFilNme'),
                                        os.getenv(key='slurmInstExt'),
                                        0o744)  # set permissions to "rwxr--r--"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--configuration", "-c", type=str, default="reticula_config.env",
                        help="Supply the full path to an environment configuration (.env) file.")
    return parser


if __name__ == '__main__':
    args = build_parser().parse_args(sys.argv[1:])
    env_path = Path(args.configuration)
    load_dotenv(dotenv_path=env_path)
    generate_program_files()
    execute_program()
