#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

guideseq.py
===========

V2

serves as the wrapper for all guideseq pipeline

"""

import os
import sys
import yaml
import argparse
import traceback
import subprocess
# Set up logger
import log
logger = log.createCustomLogger('root', "guideseq_V2")

from alignReads import alignReads
from filterBackgroundSites import filterBackgroundSites, filterBlackList, filterControl
from visualization import visualizeOfftargets
import identifyOfftargetSites
import validation
from tabulate import tabulate
from demultiplex import demultiplex, reformat_fastq

DEFAULT_YAML = os.path.dirname(os.path.realpath(__file__)) + "/default.yaml"

def get_parameters(manifest_data):

    # init
    default_refseqName = os.path.dirname(os.path.realpath(__file__)) + "/refseq_gene_name.py"
    with open(DEFAULT_YAML, 'r') as f:
        default = yaml.load(f)
    with open(manifest_data, 'r') as f:
        return_dict = yaml.load(f)  # this is user input YAML
    default['analysis_folder'] = os.getcwd()
    default['refseq_names'] = default_refseqName
    default['Manhattan.R'] = os.path.dirname(os.path.realpath(__file__)) + "/Manhattan.R"
    validation.validateManifest(return_dict)

    # assign default if user YAML missing some parameters
    return_dict = assign_default(return_dict, default)
    return_dict['control_coord_chr'], tmp = return_dict['control_primer_coord'].split(":")
    return_dict['control_coord_start'], return_dict['control_coord_end'] = [int(x) for x in tmp.split("-")]

    return return_dict

def assign_default(return_dict, default):

    section_list = ['demultiplex']
    for k in section_list:
        if k in return_dict:
            for p in default[k]:
                if not p in return_dict[k]:
                    return_dict[k][p] = default[k][p]
        else:
            return_dict[k] = default[k]

    # default global variables
    for p in default:
        if not p in section_list + ['samples']:
            if not p in return_dict:
                return_dict[p] = default[p]

    # set abs path for input and output
    for k in section_list:
        return_dict[k]['out_dir'] = os.path.join(return_dict['analysis_folder'], return_dict[k]['out_dir'])
        return_dict[k]['input_dir'] = os.path.join(return_dict['analysis_folder'], return_dict[k]['input_dir'])

    return return_dict

def get_I1_I2(R1, R2):
    """Infer I1 I2 file name

    Assumption, the only difference between R1 R2 should be just "1" or "2"

    """
    R1_list = list(R1)
    R2_list = list(R2)
    I1_list = list(R1)
    I2_list = list(R2)
    # can't directly mutate a string, have to use a list
    for i in range(len(R1_list)):
        if R1_list[i] != R2_list[i]:
            if R1_list[i - 1] == "r":
                I1_list[i - 1] = "i"
                I2_list[i - 1] = "i"
            else:
                I1_list[i - 1] = "I"
                I2_list[i - 1] = "I"
            return "".join(I1_list), "".join(I2_list)


class GuideSeq:

    def __init__(self):
        pass

    def parseManifest(self, manifest_path, sample='all', user_dict=None):
        logger.info('Loading manifest...')

        try:
            manifest_data = get_parameters(manifest_path)
            self.samples = {}
            self.parameters = manifest_data
            if user_dict == "None":
                user_dict = None
            if user_dict:
                for p in user_dict:
                    self.parameters[p] = user_dict[p]
            if sample != "all":
                self.samples[sample] = manifest_data['samples'][sample]
            else:
                self.samples = manifest_data['samples']
            del self.parameters['samples']
            logger.info("\n" + tabulate([(k, v) for k, v in self.parameters.items()]))
            logger.info("\n" + tabulate([(k, v) for k, v in self.samples[list(self.samples.keys())[0]].items()]))

        except Exception as e:
            logger.error('Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.')
            logger.error(traceback.format_exc())
            sys.exit()

    def demultiplex(self):
        barcode_dict = {}
        for sample in self.samples:
            barcode1 = self.samples[sample]['barcode1']
            barcode2 = self.samples[sample]['barcode2']
            barcode_dict[sample] = [barcode1, barcode2]

            # Check if control barcodes exist, skip if not provided
            if 'controlbarcode1' in self.samples[sample] and 'controlbarcode2' in self.samples[sample]:
                barcode1 = self.samples[sample]['controlbarcode1']
                barcode2 = self.samples[sample]['controlbarcode2']
                barcode_dict["Control_" + sample] = [barcode1, barcode2]
            else:
                logger.info(f"No control barcodes found for {sample}, skipping control processing.")

        is_file_flag = os.path.join(self.parameters['demultiplex']['input_dir'], self.parameters['demultiplex']['forward'])
        is_file_flag = os.path.isfile(is_file_flag)
        try:
            if is_file_flag:
                logger.info('Demultiplexing files...')
                demultiplex(os.path.join(self.parameters['demultiplex']['input_dir'], self.parameters['demultiplex']['forward']),
                            os.path.join(self.parameters['demultiplex']['input_dir'], self.parameters['demultiplex']['reverse']),
                            os.path.join(self.parameters['demultiplex']['input_dir'], self.parameters['demultiplex']['index1']),
                            os.path.join(self.parameters['demultiplex']['input_dir'], self.parameters['demultiplex']['index2']),
                            barcode_dict,
                            out_dir=self.parameters['demultiplex']['out_dir'],
                            mismatch=self.parameters['demultiplex']['mismatch'],
                            min_reads=self.parameters['demultiplex']['min_reads'],
                            subsample_reads=self.parameters['demultiplex']['subsample_reads'])
                logger.info('Successfully demultiplexed reads.')
            else:
                logger.error("Forward and reverse reads not found for demultiplexing.")
        except Exception as e:
            logger.error('Error demultiplexing reads.')
            logger.error(traceback.format_exc())
            quit()

    def alignReads(self):
        logger.info('Aligning reads...')
        try:
            for sample in self.samples:
                # Align reads for the treatment sample
                alignReads(self.parameters['reference_genome'],
                           self.samples[sample]['read1'],
                           self.samples[sample]['read2'],
                           os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.dedup.sam'),
                           njobs=self.parameters['njobs'],
                           bwa=self.parameters['bwa'],
                           samtools=self.parameters['samtools'],
                           umi_tools=self.parameters['umi_tools'])
                logger.info(f'Finished aligning reads for {sample}.')

                # Skip control reads if control barcodes are not provided
                if 'controlread1' in self.samples[sample] and 'controlread2' in self.samples[sample]:
                    logger.info(f'Aligning control reads for {sample}.')
                    alignReads(self.parameters['reference_genome'],
                               self.samples[sample]['controlread1'],
                               self.samples[sample]['controlread2'],
                               os.path.join(self.parameters['analysis_folder'], 'aligned', "Control_" + sample + '.dedup.sam'),
                               njobs=self.parameters['njobs'],
                               bwa=self.parameters['bwa'],
                               samtools=self.parameters['samtools'],
                               umi_tools=self.parameters['umi_tools'])
                    logger.info(f'Finished aligning control reads for {sample}.')
                else:
                    logger.info(f"No control reads for {sample}, skipping control alignment.")
        except Exception as e:
            logger.error('Error aligning reads.')
            logger.error(traceback.format_exc())
            quit()
            
    def parallel(self, manifest_path, lsf, step, overwrite):
        logger.info('Submitting parallel jobs for GuideSeqV2')
        if not os.path.exists("HPC_parallel_log"):
            os.makedirs("HPC_parallel_log")
        current_script = __file__
        count = 1
        try:
            for sample in self.samples:
                cmd = f'python {current_script} main --manifest {manifest_path} --sample {sample} --step {step} --overwrite "{overwrite}"'
                # logger.info(cmd)
                # continue
                # subprocess.call(lsf.replace("{Sample_Name}",sample).split() + [f"-J {sample[:20]}"] + [cmd])
                cmd = lsf.replace("{Sample_Name}", sample) + f" -J {sample[:10]} " + cmd
                logger.info(cmd)
                subprocess.call(cmd, shell=True)
                count += 1
            logger.info('Finished job submission')

        except Exception as e:
            logger.error('Error submitting jobs.')
            logger.error(traceback.format_exc())



def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(description='Two run modes: main or parallel',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')

    parallel_parser = subparsers.add_parser('parallel', help='submit a job for each sample specified in the YAML file. Tested in the LSF system. May not work in other job management systems.')
    parallel_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    parallel_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')
    parallel_parser.add_argument('--lsf', '-l', help='Specify LSF CMD', default='bsub -n 12 -R "span[hosts=1] rusage[mem=12000]" -P GUIDEseqV2 -q standard -o HPC_parallel_log/GUIDEseqV2_{Sample_Name}_%J.log')
    parallel_parser.add_argument('--step', help='Specify which steps of the pipeline to run (demultiplex, align, identify, visualize)', default='demultiplex+align+identify+visualize')
    parallel_parser.add_argument('--overwrite', help='overwrite specifications in the yaml file', default=None, type=yaml.load)

    main_parser = subparsers.add_parser('main', help='Run a single step or a series of steps for one or all samples')
    main_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
    main_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')
    main_parser.add_argument('--step', help='Specify steps, demultiplex, align, identify, visualize, order does not matter', default='demultiplex+align+identify+visualize')
    main_parser.add_argument('--overwrite', help='overwrite specifications in the yaml file', default=None, type=yaml.load)

    init_parser = subparsers.add_parser('init', help='initialize a default yaml in the current folder')

    return parser.parse_args()


def main():
    args = parse_args()

    if args.command == "init":
        subprocess.call(f"cp {DEFAULT_YAML} Guideseq.yaml", shell=True)
        print("Guideseq.yaml has been initialized. \nPlease add sample info and then run GUIDE-seq analysis using the 'main' or the 'parallel' subcommand.")
    elif args.command == "main":
        logger.info("User command: %s" % (" ".join(sys.argv)))
        g = GuideSeq()
        g.parseManifest(args.manifest, args.sample, args.overwrite)
        steps = args.step.split("+")
        if "demultiplex" in steps:
            g.demultiplex()
        if "align" in steps:
            g.alignReads()
        if "identify" in steps:
            g.identifyOfftargetSites()
            g.filterBackgroundSites()
        if "visualize" in steps:
            try:
                g.identified
            except:
                g.identified = {}
                for sample in g.samples:
                    if 'identified' in g.samples[sample]:
                        g.identified[sample] = g.samples[sample]['identified']
                        g.identified["Control_" + sample] = g.samples[sample]['controlidentified']
                    else:
                        g.identified[sample] = os.path.join(g.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
                        sample = "Control_" + sample
                        g.identified[sample] = os.path.join(g.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
            g.visualize()

    elif args.command == 'parallel':
        c = GuideSeq()
        c.parseManifest(args.manifest, args.sample, args.overwrite)
        # exit()
        c.parallel(args.manifest, args.lsf, args.step, args.overwrite)


if __name__ == '__main__':
    main()
