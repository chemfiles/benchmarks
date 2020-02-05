#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import gc
import sys
import tempfile
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
        self.path = "/tmp/ramdisk-" + str(os.getpid())
        os.mkdir(self.path)
        print("we need sudo access to create a ramdisk")
        check_output([
            "sudo", "mount", "-t", "tmpfs", "-o", "size={}m".format(size_in_mb), "", self.path
        ])

    def root_path(self):
        return self.path

    def close(self):
        check_output([
            "sudo", "umount", self.path
        ])
        check_output([
            "rm", "-rf", self.path
        ])


class TmpDisk:
    def __init__(self):
        self.dir = tempfile.TemporaryDirectory()

    def root_path(self):
        return self.dir.name

    def close(self):
        self.dir.cleanup()


class ram_disk:
    def __init__(self, use_ramdisk=True, size_in_mb=512):
        if use_ramdisk:
            if sys.platform.startswith("linux"):
                self.disk = LinuxRamDisk(size_in_mb)
            elif sys.platform == "darwin":
                self.disk = OsXRamDisk(size_in_mb)
            else:
                raise Exception("unable to create a RAM disk for the current OS")
        else:
            self.disk = TmpDisk()

    def __enter__(self):
        return self.disk.root_path()

    def __exit__(self, *args):
        # force gc collections to close files before atempting to close the disk
        gc.collect()
        self.disk.close()
