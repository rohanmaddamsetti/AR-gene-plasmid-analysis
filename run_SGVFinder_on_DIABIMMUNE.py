#===============================================================================
# run_SGVFinder_on_DIABIMMUNE.py by Rohan Maddamsetti.
# I modified linear_example.py bt the SGVFinder authors, for my purposes.
# Run SGVFinder on ICRA output.
# NOTE: Memory requirements are at least 20GB per file, and even then - it takes
# time, depending on the sample.
# NOTE2: See the README for requirements etc.
#===============================================================================
from glob import glob
from os.path import join, splitext, basename
from SGVFinder import get_sample_map, work_on_collection
import os
from multiprocessing import Process, Manager


## Fill in the values of the samp_to_map dict using multiprocessing.
def fill_samp_to_map(samp_to_map, f):
    ## f is the full path of the jsdel input file.
    jsdel_sample_name = splitext(basename(f))[0]
    samp_to_map[jsdel_sample_name] = get_sample_map(f, .01, 100, 10)

    
## have to run as a module in order for multiprocessing to work.
if __name__ == '__main__':

    OUTPUT_FOLDER = "/work/rm431/SGVFinder/results/DIABIMMUNE"
    try:
        os.mkdir(join(OUTPUT_FOLDER, 'browser'))
    except OSError:
        if not os.path.isdir(join(OUTPUT_FOLDER, 'browser')):
            raise

    ## I adapted the following snippet from Stack Overflow:
    ## https://stackoverflow.com/questions/38393269/fill-up-a-dictionary-in-parallel-with-multiprocessing
    manager = Manager()
    ## initialize samp_to_map here, and will update in function fill_samp_to_map
    samp_to_map = manager.dict()
    jsdel_input_files = [f for f in glob(join(OUTPUT_FOLDER, '*.jsdel'))]
    job = [Process(target=fill_samp_to_map, args=(samp_to_map, f)) for f in jsdel_input_files]
    _ = [p.start() for p in job]
    _ = [p.join() for p in job]

    ## convert the DictProxy object created by multiprocessing into a Dict.
    samp_to_map = dict(samp_to_map)
    
    ## now run the second stage of SGVFinder on the complete samp_to_map dictionary.
    vsgv, dsgv = work_on_collection(samp_to_map, 10, 0.25, 0.95, 0.25, 0.125, 0.02, 0.95, 'betaprime', 0.01, 10, 85, join(OUTPUT_FOLDER, 'browser'))
    ## and save the results to file.
    vsgv.to_pickle(join(OUTPUT_FOLDER, 'vsgv.df'))
    dsgv.to_pickle(join(OUTPUT_FOLDER, 'dsgv.df'))
    ## save as CSV files as well.
    vsgv.to_csv(join(OUTPUT_FOLDER, 'vsgv.csv'))
    dsgv.to_csv(join(OUTPUT_FOLDER, 'dsgv.csv'))
