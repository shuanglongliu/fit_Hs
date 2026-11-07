"""Set up, submit and harvest DFT jobs, one directory per spin configuration.

Each job directory of a system holds data.py (input templates and settings) and main.py (the
configurations). main.py passes the data module to the classes here and runs an action given on
the command line, e.g. ``python main.py setup``.
"""
