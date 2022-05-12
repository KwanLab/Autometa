#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Script to test autometa/common/markers.py
"""

import argparse
import pytest

from autometa import __main__


@pytest.fixture(name="no_arg_mock_parser")
def fixture_no_arg_mock_parser(monkeypatch):
    def return_mock_parser(*args, **kwargs):
        return MockParser()

    class MockParseArgs:
        citation = False
        version = False

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockParseArgs()

    # Defining the MockParser class to represent parser
    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)


@pytest.fixture(name="citation_mock_parser")
def fixture_citation_mock_parser(monkeypatch):
    def return_mock_parser(*args, **kwargs):
        return MockParser()

    class MockParseArgs:
        citation = True
        version = False

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockParseArgs()

    # Defining the MockParser class to represent parser
    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)


@pytest.fixture(name="version_mock_parser")
def fixture_version_mock_parser(monkeypatch):
    def return_mock_parser(*args, **kwargs):
        return MockParser()

    class MockParseArgs:
        citation = True
        version = False

    class MockParser:
        def add_argument(self, *args, **kwargs):
            pass

        def parse_args(self):
            return MockParseArgs()

    # Defining the MockParser class to represent parser
    monkeypatch.setattr(argparse, "ArgumentParser", return_mock_parser, raising=True)


@pytest.mark.entrypoint
def test_autometa_main_no_args(monkeypatch, no_arg_mock_parser):
    __main__.main()


@pytest.mark.entrypoint
def test_autometa_main_citation(monkeypatch, citation_mock_parser):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        __main__.main()
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0


@pytest.mark.entrypoint
def test_autometa_main_version(monkeypatch, version_mock_parser):
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        __main__.main()
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == 0
