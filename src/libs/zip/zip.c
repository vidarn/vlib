#include "miniz.h"

struct ZipFile
{
	mz_zip_archive archive;
};

struct ZipFile *zip_open_filename(const char *filename)
{
	struct ZipFile *zipfile = calloc(1, sizeof(struct ZipFile));
	mz_bool result = mz_zip_reader_init_file(&zipfile->archive, filename, 0);
	return zipfile;
}

void zip_close(struct ZipFile *zipfile)
{
	mz_bool result = mz_zip_reader_end(&zipfile->archive);
}

int zip_get_num_files(struct ZipFile *zipfile)
{
	return mz_zip_reader_get_num_files(&zipfile->archive);
}

const char *zip_get_filename(int i, struct ZipFile *zipfile)
{
	mz_zip_archive_file_stat file_stat = { 0 };
	mz_bool ret = mz_zip_reader_file_stat(&zipfile->archive, i, &file_stat);
	if (ret) {
		return file_stat.m_filename;
	}
	return 0;
}

//NOTE(Vidar):Compile miniz.c as well...

#include "miniz.c"

