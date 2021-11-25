import os

import jinja2
from jinja2 import StrictUndefined

from analyse_scan import get_path_neighbours

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


def template_fit(args, start_frame, next_frame, prev_frame=None, first=False):
    if first:
        TEMPLATE = FF_TEMPLATE
    else:
        TEMPLATE = F_TEMPLATE

    output = TEMPLATE.render(cubefit_path=args.cubefit_path, fdcm_path=args.fdcm_path, ars_path=args.ars_path,
                             n_charges=args.n_charges, suffix=args.suffix,
                             n_steps=args.n_steps, scan_name=args.scan_name, cubes_dir=args.cubes_dir,
                             output_dir=args.output_dir, frames=args.frames, initial_fit=args.initial_fit,
                             initial_fit_cube=args.initial_fit_cube, n_scan_points=args.n_scan_points,
                             start_frame=start_frame, next_frame=next_frame, prev_frame=prev_frame, acd=args.acd)
    return output


def template_neighbours(args, do_neighbours=True):
    paths, neighbours = get_path_neighbours(args)

    if not os.path.exists(args.job_folder):
        os.makedirs(args.job_folder)

    n_jobs = len(paths)

    for i, (path, neighbour) in enumerate(zip(paths, neighbours)):
        if i < n_jobs - 1:
            is_first = (i == 0)
            print(i, path, neighbour)
            tmp_str = template_fit(args, paths[i], paths[i + 1], first=is_first)
            f = open(os.path.join(args.job_folder, f"frame_{paths[i]}_{paths[i + 1]}.sh"), "w")
            f.write(tmp_str)
            if do_neighbours:
                for n in neighbour:
                    tmp_str = template_fit(args, paths[i], n)
                    _fpath = os.path.join(args.job_folder, f"frame_{paths[i]}_{n}.sh")
                    f_ = open(_fpath, "w")
                    f_.write(tmp_str)
                    f_.close()
                    f.write(f"\nsbatch {_fpath} \n")

            try:
                next_job = os.path.join(args.job_folder, f"frame_{paths[i + 1]}_{paths[i + 2]}.sh")
                f.write(f"\nsbatch {next_job} \n")
            except IndexError:
                pass

            f.close()
