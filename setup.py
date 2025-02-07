#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import guideseq

import os
if os.path.isfile("README.md"):
	with open("README.md", "r") as fh:
		long_description = fh.read()
else:
	long_description="guide-seq"


setup(
	name='guide_seq',
	version=str(guideseq.__version__),
	description="GUIDE-seq pipeline for GUIDE-seq and GUIDE-seq-2",
	author="Yichao Li, Shengdar Q Tsai",
	author_email='Yichao.Li@stjude.org',
	url='https://github.com/tsailabSJ/guideseq',
	# packages=find_packages(),
	packages=[
		'guideseq',
	],
	package_dir={'guideseq':
				 'guideseq'},
	
	scripts=['guideseq/guideseq.py','guideseq/alignReads.py','guideseq/visualization.py',
		'guideseq/filterBackgroundSites.py','guideseq/identifyOfftargetSites.py'],
	package_data={'test': ["test/*"]},
	license="AGPL",
	include_package_data=True,
	long_description=long_description,
	long_description_content_type='text/markdown',
	keywords='guideseq',
	classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'Topic :: Scientific/Engineering :: Visualization',
		'Topic :: Scientific/Engineering :: Information Analysis',
		'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
		'Operating System :: Unix',
		'Natural Language :: English',
		'Programming Language :: Python :: 3'
	]
)
