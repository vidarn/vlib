#include "data_buffer.h"
#include <stdlib.h>
#include <string.h>

//TODO(Vidar): Do we want versions for each entry??

unsigned char * data_buffer_to_stream(struct DataBuffer data, int *size)
{
    *size = data.total_size + sizeof(struct DataBuffer);
    unsigned char *ret = calloc(*size,1);
    struct DataBuffer *tmp = (struct DataBuffer*)ret;
    *tmp = data;
    memcpy(tmp+1, data.data, data.total_size);
    return ret;
}

struct DataBuffer data_buffer_from_stream(unsigned char *data)
{
    struct DataBuffer *d = (struct DataBuffer *)data;
    struct DataBuffer ret = *d;
    ret.data = calloc(ret.total_size, 1);
    memcpy(ret.data,d+1,ret.total_size);
    return ret;
}


/* Entry layout:
   int size,
   int name_size (including null)
   char name[name_size]
   unsigned char data[size]
 */

void data_buffer_add_value(const char *name, void *val, int size,
    struct DataBuffer *data)
{
    int str_size = (int)strlen(name) + 1;
    int old_size = data->total_size;
    int new_size = data->total_size + size + str_size + sizeof(int) * 2;
    data->data = realloc(data->data, new_size);
    int *entry_int = (int*)(data->data+old_size);
    entry_int[0] = size;
    entry_int[1] = str_size;
    char * entry_char = (char *)(entry_int+2);
    memcpy(entry_char, name, str_size);
    entry_char += str_size;
    memcpy(entry_char, val, size);
    data->num_entries++;
    data->total_size = new_size;
}

struct DataBufferDecode data_buffer_decode_begin(struct DataBuffer data)
{
    struct DataBufferDecode ret = {0};
    ret.buffer = data;
    return ret;
}

struct DataBufferEntry data_buffer_decode_get_value(struct DataBufferDecode *d)
{
    struct DataBufferEntry entry = {0};
    if(d->buffer.num_entries > 0){
        int *tmp_int = (int*)d->buffer.data;
        int size     = tmp_int[0];
        int str_size = tmp_int[1];
        int entry_size = size + str_size + 2*sizeof(int);
        unsigned char *tmp_char = (unsigned char *)(tmp_int + 2);
        entry.name = (char*)tmp_char;
        tmp_char += str_size;
        entry.data = tmp_char;
        d->buffer.num_entries--;
        d->buffer.total_size -= entry_size;
        d->buffer.data += entry_size;
        entry.type = DATA_BUFFER_ENTRY_TYPE_VALUE;
    }else{
        entry.type = DATA_BUFFER_ENTRY_TYPE_EOF;
    }
    return entry;
}
