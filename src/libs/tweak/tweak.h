#pragma once

enum TweakValueType{
    TWEAK_VALUE_TYPE_FLOAT,
    TWEAK_VALUE_TYPE_INT,
};

struct TweakValue{
    char *name;
    union{
        float val_float;
        int val_int;
    };
    enum TweakValueType type;
};

struct TweakValueCollection{
    int num_tweak_values;
    struct TweakValue *values;
};

float tweak_float(const char *name, float default_value,
    struct TweakValueCollection *collection);
int tweak_int(const char *name, int default_value,
                  struct TweakValueCollection *collection);

#include "gui/gui.h"
void tweak_values_gui(struct TweakValueCollection *collection,
    struct GuiContext *gui);
