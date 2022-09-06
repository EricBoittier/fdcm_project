import os
import pandas as pd
import jinja2
from jinja2 import StrictUndefined

from create_neighbours import schedule_from_2drmsd, scan_neighbours_schedule_from_anytree_graph

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
                               morton=args.morton, morton_start=args.morton_start, acd=args.acd)
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


def template_neighbours_from_2drmsd(args, do_neighbours=True):
    scan_neighbours_schedule = schedule_from_2drmsd(args.rmsd)
    #  Make directory if needed
    if not os.path.exists(args.job_folder):
        os.makedirs(args.job_folder)

    #  list to keep track of jobs so a file isn't overwritten.
    jobs_submitted = []
    #  Process the schedule for job writing
    for i, (scan, neighbours, schedule) in enumerate(scan_neighbours_schedule):
        # print(i, scan, neighbours, schedule)

        previous, start, end = scan
        is_first = (start == previous)
        tmp_str = template_fit(args, start, end, first=is_first, prev_frame=previous)
        with open(os.path.join(args.job_folder, f"frame_{start}_{end}.sh"), "w") as f:
            jobs_submitted.append(os.path.join(args.job_folder, f"frame_{start}_{end}.sh"))
            print(os.path.join(args.job_folder, f"frame_{start}_{end}.sh"))
            f.write(tmp_str)

            #  write the neighbour jobs and schedule them
            if do_neighbours:
                for n in neighbours:
                    _fpath = os.path.join(args.job_folder, f"frame_{start}_{n[1]}.sh")
                    if _fpath not in jobs_submitted:
                        tmp_str = template_fit(args, start, n[1], first=is_first, prev_frame=previous)

                        f_ = open(_fpath, "w")
                        f_.write(tmp_str)

                        print(f"\nsbatch {_fpath} \n")
                        f.write(f"\nsbatch {_fpath} \n")

            #  schedule the jobs preceding this one (i.e. all branches from this node)
            for job in schedule:
                next_job = os.path.join(args.job_folder, f"frame_{job[0]}_{job[1]}.sh")
                f.write(f"\nsbatch {next_job} \n")
                print(next_job)


def template_neighbours_from_anytree_and_G(args, do_neighbours=False):
    anytree = pd.read_pickle(args.anytree)
    print(anytree)
    G = pd.read_pickle(args.g_obj)

    scan_neighbours_schedule = scan_neighbours_schedule_from_anytree_graph(anytree, G)

    #  Make directory if needed
    if not os.path.exists(args.job_folder):
        os.makedirs(args.job_folder)

    #  list to keep track of jobs so a file isn't overwritten.
    jobs_submitted = []
    #  Process the schedule for job writing
    for i, (scan, neighbours, schedule) in enumerate(scan_neighbours_schedule):
        print(i, scan, neighbours, schedule)

        previous, start, end = scan
        is_first = (0 == previous)
        tmp_str = template_fit(args, start, end, first=is_first, prev_frame=previous)
        with open(os.path.join(args.job_folder, f"p{start}_{end}.sh"), "w") as f:
            jobs_submitted.append(os.path.join(args.job_folder, f"p{start}_{end}.sh"))
            print(os.path.join(args.job_folder, f"p{start}_{end}.sh"))
            f.write(tmp_str)

            #  write the neighbour jobs and schedule them
            if do_neighbours:
                for n in neighbours:
                    _fpath = os.path.join(args.job_folder, f"p{start}_{n[1]}.sh")
                    if _fpath not in jobs_submitted:
                        tmp_str = template_fit(args, start, n[1], first=is_first, prev_frame=previous)

                        f_ = open(_fpath, "w")
                        f_.write(tmp_str)

                        print(f"\nsbatch {_fpath} \n")
                        f.write(f"\nsbatch {_fpath} \n")

            #  schedule the jobs preceding this one (i.e. all branches from this node)
            for job in schedule:
                next_job = os.path.join(args.job_folder, f"p{job[0]}_{job[1]}.sh")
                f.write(f"\nsbatch {next_job} \n")
                print(next_job)
