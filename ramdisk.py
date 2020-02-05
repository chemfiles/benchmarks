#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
from subprocess import check_output


class OsXRamDisk:
    def __init__(self, size_in_mb):
        blocksize = 2048 * size_in_mb
        output = check_output([
            "hdiutil", "attach", "-nomount", "ram://{}".format(blocksize)
        ])
        self.hdid = output.strip()
        self.path = "ramdisk-" + str(os.getpid())
        check_output([
            "diskutil", "partitionDisk", self.hdid, "1", "GPTFormat", "APFS", self.path, "100%"
        ])

    def root_path(self):
        return "/Volumes/{}".format(self.path)

    def close(self):
        check_output([
            "diskutil", "eject", self.hdid
        ])


class LinuxRamDisk:
    def __init__(self, size_in_mb):
        pass

    def close(self):
        pass


class ram_disk:
    def __init__(self, size_in_mb=512):
        if sys.platform == "linux2":
            self.disk = LinuxRamDisk(size_in_mb)
        elif sys.platform == "darwin":
            self.disk = OsXRamDisk(size_in_mb)
        else:
            raise Error("unable to create a RAM disk for the current OS")

    def __enter__(self):
        return self.disk.root_path()

    def __exit__(self, *args):
        self.disk.close()
