import os, fnmatch

def get_files_in_traverse_dir(a_dir, a_file):
    """ return list of all files in directory and its subdirectories \
    which matches 'a_file' and its subdirectory path, support Unix \
    filename pattern matching ('*', '?', [seq], [!seq]) """

    ls_file = []
    for root, dirs, files in os.walk(a_dir):
        for name in files:
            if fnmatch.fnmatch(name, a_file):
                subdir = root.replace(a_dir, '')
                ls_file.append((os.path.join(root, name), subdir))
    return ls_file


def sort_get_z(files, a_sim_info):
    zs = []
    for a_file in files:
        if a_sim_info.app + '_z' in a_file:
            z = float(a_file[a_file.index(a_sim_info.app + '_z') + len(a_sim_info.app+'_z'):-4])
        elif a_sim_info.app + '_init' in a_file:
            z = 'init'
        else:
            print "WARNING! Skipping file '%s', unknown format." % a_file
        zs.append(z)
    return zip(*sorted(zip(zs, files), reverse=True))

def sort_get_fl_get_z(a_sim_info, subdir, a_file='*.dat'):
    files = [x[0] for x in get_files_in_traverse_dir(a_sim_info.dir + subdir, a_file)]
    return sort_get_z(files, a_sim_info)

def slice_zs_files(zs, files, a_slice=2.4):
    a_ = 0
    files_sliced = []
    zs_sliced = []
    for i, z in enumerate(zs):
        a = 1/(1+z)
        if (a < a_slice*a_) and a != 1: continue
        a_ = a
        zs_sliced.append(z)
        files_sliced.append(files[i])
    return zs_sliced, files_sliced

def create_dir(out_dir):
    if not os.path.exists(out_dir):
        # print "Creating outdir '%s'" % out_dir
        os.makedirs(out_dir)

def try_get_zs_files(a_sim_info, subdir, a_file='*.dat'):
    try:
        zs, files = sort_get_fl_get_z(a_sim_info, subdir, a_file=a_file)
    except ValueError:
        print "WARNING! Missing data in '%s'. Skipping step." % (a_sim_info.dir + subdir)
        return None, None
    else:
        return list(zs), list(files)
