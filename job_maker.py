import os

import jinja2
from jinja2 import StrictUndefined

from analyse_scan import get_path_neighbours_from_gaussian_scan
from create_neighbours import schedule_from_2drmsd

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


def template_neighbours_from_gaussian_scan(args, do_neighbours=True):
    paths, neighbours = get_path_neighbours_from_gaussian_scan(args)
    if not os.path.exists(args.job_folder):
        os.makedirs(args.job_folder)
    n_jobs = len(paths)
    for i, (path, neighbour) in enumerate(zip(paths, neighbours)):
        if i < n_jobs - 1:
            is_first = (i == 0)
            print(i, path, neighbour)
            tmp_str = template_fit(args, paths[i], paths[i + 1], first=is_first, prev_frame=paths[i - 1])
            f = open(os.path.join(args.job_folder, f"frame_{paths[i]}_{paths[i + 1]}.sh"), "w")
            f.write(tmp_str)
            if do_neighbours:
                for n in neighbour:
                    tmp_str = template_fit(args, paths[i], n, first=is_first, prev_frame=paths[i - 1])
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


def template_neighbours_from_2drmsd(args, do_neighbours=True):
    scan_neighbours_schedule = schedule_from_2drmsd(args.rmsd)
    #  Make directory if needed
    if not os.path.exists(args.job_folder):
        os.makedirs(args.job_folder)
    #  Process the schedule for job writing
    for i, (scan, neighbours, schedule) in enumerate(scan_neighbours_schedule):
        print(i, scan, neighbours, schedule)
        is_first = (i == 0)
        previous, start, end = scan
        tmp_str = template_fit(args, start, end, first=is_first, prev_frame=previous)

        f = open(os.path.join(args.job_folder, f"frame_{start}_{end}.sh"), "w")
        f.write(tmp_str)
        #  write the neighbour jobs and schedule them
        if do_neighbours:
            for n in neighbours:
                tmp_str = template_fit(args, start, n[1], first=is_first, prev_frame=previous)
                _fpath = os.path.join(args.job_folder, f"frame_{start}_{n[1]}.sh")
                f_ = open(_fpath, "w")
                f_.write(tmp_str)
                f_.close()
                print(f"\nsbatch {_fpath} \n")
                f.write(f"\nsbatch {_fpath} \n")

        #  schedule the jobs preceding this one (i.e. all branches from this node)
        for job in schedule:
            next_job = os.path.join(args.job_folder, f"frame_{job[0]}_{job[1]}.sh")
            f.write(f"\nsbatch {next_job} \n")
            print(next_job)

        f.close()
