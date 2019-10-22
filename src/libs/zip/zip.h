#pragma once

struct ZipFile
;

struct ZipFile *zip_open_filename(const char *filename)
;

void zip_close(struct ZipFile *zipfile)
;

int zip_get_num_files(struct ZipFile *zipfile)
;

const char *zip_get_filename(int i, struct ZipFile *zipfile)
;