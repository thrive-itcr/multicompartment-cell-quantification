# Copyright (c) General Electric Company, 2017.  All rights reserved.

# Rt 106

# Multi-Compartment Cell Quantification

import os, glob, uuid, time, logging, subprocess

# function: run_algorithm() -- Python function for marshalling your data and running your algorithm.
# parameters:
#   datastore: object tobe used when interacting with the Object Store
#   context:  A JSON structure that should contain all the inputs and parameters your algorithm needs.
def run_algorithm(datastore,context):

    logging.info('run_algorithm: %r' % context)

    # cleanup the input and output directories
    for f in glob.glob('/rt106/input/*') + glob.glob('/rt106/output/*'):
        os.remove(f)

# 1.    Code for marshalling inputs.
    nuc_seg_path =  datastore.get_pathology_result_image_path(context['slide'], context['region'], context['branch'], 'NucSeg')
    nuc_seg_image = 'NucSeg.tif'
    instance_status = datastore.get_instance(nuc_seg_path, '/rt106/input', nuc_seg_image, 'tiff16')
    # Check for error.
    if (instance_status != 200 ):
        status = "ERROR_NucSeg_FILE_NOT_FOUND"
        return { 'result' : {}, 'status' : status }

    
    if 'quantitateCells' in context and context['quantitateCells']:
        cell_seg_path = datastore.get_pathology_result_image_path(context['slide'], context['region'], context['branch'],'CellSeg')
        cell_seg_image = 'CellSeg.tif'
        instance_status = datastore.get_instance(cell_seg_path, '/rt106/input', cell_seg_image, 'tiff16')
        # Check for error.  The line below tests for an empty list of instances having been returned.
        if instance_status != 200:
            status = "ERROR_CellSeg_FILE_NOT_FOUND"
            return { 'result' : {}, 'status' : status }

    datastore.retrieve_multi_channel_pathology_image(context['slide'], context['region'], '/rt106/input')

    output_path = datastore.get_pathology_result_path(context['slide'], context['region'], context['branch'], 'Quant')
    logging.info('output_path: %r' % output_path)
    quan_csv = 'quant_%s.csv' % context['region']
    output_file = '/rt106/output/%s' % quan_csv

# 2.    Code for calling algorithm.
    try:
		#Parameters for multiple channel quantification
		# usage: -insegmask SegmentationMask -inbioim BioMrkerImagenames (space separated) -bioname BioNames (space separated) -outname OutFName
        bio_images = ''
        bio_names = ''
        for f in next(os.walk('/rt106/input/'))[2]:
            logging.info('filename in /rt106/input/: %s' % f)
            if f != 'NucSeg.tif' and f != 'CellSeg.tif':
                bio_images = bio_images + '/rt106/input/' + f + ' '
                logging.info('bio_images %s' % bio_images)
                bio_names = bio_names + f.split('.')[0] + ' '
                logging.info('bio_names %s' % bio_names)

        output_args = '-outname %s' % output_file
        if 'quantitateCells' in context and context['quantitateCells']:
            input_args = '-Nucsegmask /rt106/input/NucSeg.tif -inbioim %s -bioname %s %s -Cellsegmask /rt106/input/CellSeg.tif ' % (bio_images, bio_names, output_args)
        else:
            input_args = '-Nucsegmask /rt106/input/NucSeg.tif -inbioim %s -bioname %s %s ' % (bio_images, bio_names, output_args)

        logging.info('input_args %s' % input_args)
        run_algorithm = '/usr/bin/python CellQuantificationMultiMarkers.py %s' % (input_args)
        logging.info('run Algorithm: %r' % run_algorithm)
        subprocess.check_call(run_algorithm,shell=True)
    except subprocess.CalledProcessError, e:
        logging.error('%d - %s' % (e.returncode, e.cmd))
        status = "EXECUTION_FINISHED_ERROR"
        result_context = {}
        return { 'result' : result_context, 'status' : status }

# 3.    Set status.
    status = "EXECUTION_FINISHED_SUCCESS"

# 4.    Store results in datastore.
    response_upload = datastore.post_instance(output_path, '/rt106/output', quan_csv, 'csv', context['force'])
    
    if response_upload == 403:
        status = "EXECUTION_ERROR"

# 5.    Create JSON structure containing results.
    nuclear_image_path = datastore.get_pathology_primary_path(context['slide'], context['region'], 'DAPI')
    result_context = {
        "nuclearImage" : nuclear_image_path,
        "cellMetrics" : output_path
    }

    return { 'result' : result_context, 'status' : status }
