#include "tweak.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static struct TweakValue *get_tweak_val(const char *name, enum TweakValueType type,
    void *default_value, struct TweakValueCollection *collection)
{
    for(int i_value = 0; i_value < collection->num_tweak_values; i_value++){
        struct TweakValue *val = collection->values + i_value;
        if(strcmp(name, val->name)==0){
            return val;
        }
    }
    //NOTE(Vidar): If we get here, the value was not found..
    int i = collection->num_tweak_values;
    collection->num_tweak_values++;
    collection->values = realloc(collection->values, collection->num_tweak_values*sizeof(struct TweakValue));
    struct TweakValue *val= collection->values+i;
    val->name = strdup(name);
    val->type = type;
    switch (type) {
        case TWEAK_VALUE_TYPE_FLOAT:
            val->val_float = *(float*)default_value;
            break;
        case TWEAK_VALUE_TYPE_INT:
            val->val_int = *(int*)default_value;
            break;
    }
    return val;
}

float tweak_float(const char *name, float default_value,
                  struct TweakValueCollection *collection)
{
    return get_tweak_val(name, TWEAK_VALUE_TYPE_FLOAT, &default_value, collection)->val_float;
}

int tweak_int(const char *name, int default_value,
                  struct TweakValueCollection *collection)
{
    return get_tweak_val(name, TWEAK_VALUE_TYPE_INT, &default_value, collection)->val_int;
}

void tweak_values_gui(struct TweakValueCollection *collection,
    struct GuiContext *gui)
{
    struct nk_context *ctx = &gui->ctx;
    if(nk_tree_push(ctx, NK_TREE_TAB, "Tweak values", NK_MAXIMIZED)){
        nk_layout_row_dynamic(ctx, 20, 1);
        for(int i_value = 0; i_value < collection->num_tweak_values; i_value++){
            struct TweakValue *val = collection->values + i_value;
            switch (val->type) {
                case TWEAK_VALUE_TYPE_FLOAT:
                {
                    nk_property_float(ctx, val->name, 0.f, &val->val_float, 1.f, 0.1f, 0.01f);
                    break;
                }
                case TWEAK_VALUE_TYPE_INT:
                {
                    nk_property_int(ctx, val->name, 0, &val->val_int, 100000, 1.f, 1.f);
                    break;
                }
                default:
                    break;
            }
        }
        nk_tree_pop(ctx);
    }
}
