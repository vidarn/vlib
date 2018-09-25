#pragma once

struct DataBuffer{
    int num_entries;
    int total_size;
    unsigned char *data;
};

struct DataBufferDecode{
    struct DataBuffer buffer;
};

enum DataBufferEntryType {
    DATA_BUFFER_ENTRY_TYPE_NONE,
    DATA_BUFFER_ENTRY_TYPE_EOF,
    DATA_BUFFER_ENTRY_TYPE_VALUE,
};

struct DataBufferEntry{
    enum DataBufferEntryType type;
    char *name;
    int size;
    void *data;
};

unsigned char * data_buffer_to_stream(struct DataBuffer data, int *size);
struct DataBuffer data_buffer_from_stream(unsigned char *data);
void data_buffer_add_value(const char *name, void *val, int size,
    struct DataBuffer *data);
struct DataBufferDecode data_buffer_decode_begin(struct DataBuffer data);
struct DataBufferEntry data_buffer_decode_get_value(struct DataBufferDecode *d);
