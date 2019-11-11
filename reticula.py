import argparse
import os
import sys
from string import Template


def execute_program(configuration_kv_map):
    # run the generated build script
    pass


def compose_template(template_path, base_path, base_permissions, configuration_kv_map):
    with open(template_path, 'r') as build_template_filepointer:
        build_template = build_template_filepointer.read()
        build_script = Template(build_template).substitute(configuration_kv_map)
    with open(base_path, 'w') as build_script_fp:
        build_script_fp.write(build_script)
        os.chmod(base_path, base_permissions)


def generate_program_files(configuration_kv_map):
    build_tmpl = f'{configuration_kv_map.TEMPL_DIR}/{configuration_kv_map.BUILD_FN}.{configuration_kv_map.TEMPL_EXT}'
    build_base = f'{configuration_kv_map.BASE_DIR}/{configuration_kv_map.BUILD_FN}.{configuration_kv_map.SCRPT_EXT}'

    dockerfile_tmpl = f'{configuration_kv_map.TEMPL_DIR}/{configuration_kv_map.DOCKR_FN}.{configuration_kv_map.TEMPL_EXT}'
    dockerfile_base = f'{configuration_kv_map.TEMPL_DIR}/Dockerfile'  # the default dockerfile filename

    slurm_tmpl = f'{configuration_kv_map.TEMPL_DIR}/{configuration_kv_map.SLURM_FN}.{configuration_kv_map.TEMPL_EXT}'
    slurm_base = f'{configuration_kv_map.TEMPL_DIR}/{configuration_kv_map.SLURM_FN}.{configuration_kv_map.SCRPT_EXT}'

    compose_template(build_tmpl, build_base, 0o555, configuration_kv_map)  # set permissions to "r-xr-xr-x"
    compose_template(dockerfile_tmpl, dockerfile_base, 0o444, configuration_kv_map)  # set permissions to "r--r--r--"
    compose_template(slurm_tmpl, slurm_base, 0o555, configuration_kv_map)  # set permissions to "r-xr-xr-x"


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--configuration", "-c", type=str, default="reticula_config.env",
                        help="Supply the full path to a reticula_config.env configuration file.")
    return parser


def main():
    args = build_parser().parse_args(sys.argv[1:])
    configuration_kv_map = json.load(args.configuration)
    generate_program_files(configuration_kv_map)
    execute_program(configuration_kv_map)
