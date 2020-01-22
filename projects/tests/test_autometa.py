#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest
import shutil


class AutometaTestCase(unittest.TestCase):
    def setUp(self):
        self.metagenome = Metagenome('The Metagenome')

    def tearDown(self):
        shutil.rmtree(self.outdir)

    # TODO: Finish AutometaTests Class
