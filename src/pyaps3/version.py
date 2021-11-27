#!/usr/bin/env python3
# Author: Zhang Yunjun, Nov 2021
# Copyright 2012, by the California Institute of Technology.


import collections

# release history
Tag = collections.namedtuple('Tag', 'version date')
release_history = (
    Tag('0.3.1', '2021-11-26'),
    Tag('0.3.0', '2021-11-15'),
    Tag('0.2.0', '2021-08-31'),
    Tag('0.1.0', '2019-04-07'),
)

# latest release
release_version = 'v{}'.format(release_history[0].version)
release_date = release_history[0].date

