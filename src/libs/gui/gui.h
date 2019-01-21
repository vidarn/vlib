#pragma once
#define NK_INCLUDE_DEFAULT_ALLOCATOR
#define NK_INCLUDE_STANDARD_VARARGS
#define NK_INCLUDE_STANDARD_IO
#include "nuklear.h"
#include "engine.h"

struct GuiContext
{
    struct nk_context ctx;
    struct GameData *data;
    
    int rect_shader;
    int circle_shader;
    int triangle_shader;
    int color_rect_shader;
};

struct GuiContext *gui_init(int font, struct GameData *data);
void gui_begin_frame(struct GuiContext *gui, struct InputState input_state, struct GameData *data);
void gui_draw(struct GuiContext *gui, struct RenderContext *context);


void gui_fps_section(struct GuiContext *gui, struct GameData *data);
struct Color gui_color_picker(struct Color col, struct GuiContext *gui,
    struct GameData *data);
