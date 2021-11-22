import os

import jinja2
from jinja2 import StrictUndefined

"""
Templates
"""
C_TEMPLATE_STR = open(os.path.join("job_templates", "TEMPLATE_concerted.sh")).read()
C_TEMPLATE = jinja2.Template(C_TEMPLATE_STR, undefined=StrictUndefined)

S_TEMPLATE_STR = open(os.path.join("job_templates", "TEMPLATE_scan.sh")).read()
S_TEMPLATE = jinja2.Template(S_TEMPLATE_STR, undefined=StrictUndefined)

M_TEMPLATE_STR = open(os.path.join("job_templates", "TEMPLATE_morton.sh")).read()
M_TEMPLATE = jinja2.Template(M_TEMPLATE_STR, undefined=StrictUndefined)

FF_TEMPLATE_STR = open(os.path.join("job_templates", "TEMPLATE_first_fit.sh")).read()
FF_TEMPLATE = jinja2.Template(FF_TEMPLATE_STR, undefined=StrictUndefined)

F_TEMPLATE_STR = open(os.path.join("job_templates", "TEMPLATE_following_fit.sh")).read()
F_TEMPLATE = jinja2.Template(F_TEMPLATE_STR, undefined=StrictUndefined)


def template_scan(args):
    output = C_TEMPLATE.render(cubefit_path=args.cubefit_path, fdcm_path=args.fdcm_path, ars_path=args.ars_path,
                               n_charges=args.n_charges, suffix=args.suffix,
                               n_steps=args.n_steps, scan_name=args.scan_name, cubes_dir=args.cubes_dir,
                               output_dir=args.output_dir, frames=args.frames, initial_fit=args.initial_fit,
                               initial_fit_cube=args.initial_fit_cube)
    return output


def template_concerted(args):
    output = S_TEMPLATE.render(cubefit_path=args.cubefit_path, fdcm_path=args.fdcm_path, ars_path=args.ars_path,
                               n_charges=args.n_charges, suffix=args.suffix,
                               n_steps=args.n_steps, scan_name=args.scan_name, cubes_dir=args.cubes_dir,
                               output_dir=args.output_dir, frames=args.frames, initial_fit=args.initial_fit,
                               initial_fit_cube=args.initial_fit_cube, n_scan_points=args.n_scan_points)
    return output


def template_morton(args):
    output = M_TEMPLATE.render(cubefit_path=args.cubefit_path, fdcm_path=args.fdcm_path, ars_path=args.ars_path,
                               n_charges=args.n_charges, suffix=args.suffix,
                               n_steps=args.n_steps, scan_name=args.scan_name, cubes_dir=args.cubes_dir,
                               output_dir=args.output_dir, frames=args.frames, initial_fit=args.initial_fit,
                               initial_fit_cube=args.initial_fit_cube, n_scan_points=args.n_scan_points,
                               morton=args.morton, morton_start=args.morton_start)
    return output


def template_first_fit(args):
    output = FF_TEMPLATE.render(cubefit_path=args.cubefit_path, fdcm_path=args.fdcm_path, ars_path=args.ars_path,
                                n_charges=args.n_charges, suffix=args.suffix,
                                n_steps=args.n_steps, scan_name=args.scan_name, cubes_dir=args.cubes_dir,
                                output_dir=args.output_dir, frames=args.frames, initial_fit=args.initial_fit,
                                initial_fit_cube=args.initial_fit_cube, n_scan_points=args.n_scan_points,
                                start_frame=args.start_frame, next_frame=args.next_frame)
    return output


def template_following_fit(args):
    output = F_TEMPLATE.render(cubefit_path=args.cubefit_path, fdcm_path=args.fdcm_path, ars_path=args.ars_path,
                               n_charges=args.n_charges, suffix=args.suffix,
                               n_steps=args.n_steps, scan_name=args.scan_name, cubes_dir=args.cubes_dir,
                               output_dir=args.output_dir, frames=args.frames, initial_fit=args.initial_fit,
                               initial_fit_cube=args.initial_fit_cube, n_scan_points=args.n_scan_points,
                               start_frame=args.start_frame, next_frame=args.next_frame)
    return output
