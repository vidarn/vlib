#define NK_IMPLEMENTATION
#include "gui.h"
#include <stdio.h>
#include "os/win32/key_codes.h"

//TODO(Vidar):This is only for the definition of SATIN_SHADER_SUFFIX
#include "vn_gl.h"

struct GuiFont {
    int font;
    struct GuiContext *context;
};

float gui_text_width(nk_handle handle, float height, const char *text, int len)
{
    struct GuiFont *gui_font = handle.ptr;
    float text_width = get_string_render_width(gui_font->font, text, len,
        gui_font->context->data);
    return text_width*(float)reference_resolution;
}

struct GuiContext *gui_init(int font, struct GameData *data)
{
    struct GuiContext *context = calloc(1, sizeof(struct GuiContext));
    struct nk_context *ctx = &context->ctx;
    struct nk_user_font *nk_font = calloc(1, sizeof(struct nk_user_font));
    struct GuiFont *gui_font = calloc(1, sizeof(struct GuiFont));
    context->data = data;
    gui_font->font = font;
    gui_font->context = context;
    nk_font->userdata.ptr = gui_font;
    nk_font->height = get_font_height(font, data);
    nk_font->width = gui_text_width;
    nk_init_default(ctx, nk_font);
    
    context->rect_shader = load_shader("gui" SATIN_SHADER_SUFFIX, "gui_rect" SATIN_SHADER_SUFFIX, data);
    context->circle_shader = load_shader("gui" SATIN_SHADER_SUFFIX, "gui_circle" SATIN_SHADER_SUFFIX, data);
    context->triangle_shader = load_shader("gui" SATIN_SHADER_SUFFIX, "gui_triangle" SATIN_SHADER_SUFFIX, data);
    context->color_rect_shader = load_shader("gui" SATIN_SHADER_SUFFIX, "gui_color_rect" SATIN_SHADER_SUFFIX, data);
    
    return context;
}

void gui_begin_frame(struct GuiContext *gui, struct InputState input_state, struct GameData *data)
{
    struct nk_context *ctx = &gui->ctx;

    float mx = input_state.mouse_x*(float)reference_resolution;
    float my = (1.f-input_state.mouse_y)*(float)reference_resolution;
	float dx = input_state.delta_x*(float)reference_resolution;
	float dy = input_state.delta_y*(float)reference_resolution;

    nk_input_begin(ctx);
    //nk_input_motion(ctx, mx, my);
    ctx->input.mouse.pos.x = mx;
    ctx->input.mouse.pos.y = my;
    ctx->input.mouse.delta.x = dx;
	ctx->input.mouse.delta.y = dy;
    ctx->input.mouse.prev.x = mx-dx;
    ctx->input.mouse.prev.y = my-dy;
    nk_input_button(ctx, NK_BUTTON_LEFT, mx, my, input_state.mouse_down);
    //TODO(Vidar): Provide more input info

	for (int i = 0; i < input_state.num_keys_typed; i++) {
		uint32_t key = input_state.keys_typed[i];
		nk_input_unicode(ctx, key);
	}
	nk_input_key(ctx, NK_KEY_ENTER, is_key_down(KEY_RETURN));
	nk_input_key(ctx, NK_KEY_BACKSPACE, is_key_down(KEY_BACKSPACE));
	nk_input_key(ctx, NK_KEY_LEFT, is_key_down(KEY_LEFT));
	nk_input_key(ctx, NK_KEY_RIGHT, is_key_down(KEY_RIGHT));
	nk_input_key(ctx, NK_KEY_UP, is_key_down(KEY_UP));
	nk_input_key(ctx, NK_KEY_DOWN, is_key_down(KEY_DOWN));
	nk_input_key(ctx, NK_KEY_DEL, is_key_down(KEY_DELETE));
	nk_input_key(ctx, NK_KEY_TAB, is_key_down(KEY_TAB));
	nk_input_key(ctx, NK_KEY_TEXT_END, is_key_down(KEY_ESCAPE));

    nk_input_end(&gui->ctx);

	if (ctx->input.mouse.grabbed) {
		lock_cursor(data);
	}
	else {
		unlock_cursor(data);
	}

}

static float get_coord(short c)
{
    const float inv_reference_resolution = 1.f/(float)reference_resolution;
    return (float)c*inv_reference_resolution;
}
static struct Color get_color(struct nk_color c)
{
    struct Color color = {(float)c.r/255.f, (float)c.g/255.f, (float)c.b/255.f,
        (float)c.a/255.f};
    return color;
}
void gui_draw(struct GuiContext *gui, struct RenderContext *context)
{
    const struct nk_command *cmd = 0;
    float scissor_x1 = 0.f;
    float scissor_y1 = 0.f;
    float scissor_x2 = 1.f;
    float scissor_y2 = 1.f;
    
    nk_foreach(cmd, &gui->ctx) {
#define NK_COMMAND_CASE_BEGIN(NAME, name) case NK_COMMAND_##NAME:{\
struct nk_command_##name *name = (struct nk_command_##name *)cmd;

#define NK_COMMAND_CASE_END break;}

#define SCISSOR(val, dim) {\
    if(val < scissor_##dim##1) val = scissor_##dim##1;\
    if(val > scissor_##dim##2) val = scissor_##dim##2;\
}
        switch (cmd->type) {
        NK_COMMAND_CASE_BEGIN(SCISSOR, scissor)
            scissor_x1 = get_coord(scissor->x);
            scissor_y1 = get_coord(scissor->y);
            scissor_x2 = scissor_x1 + get_coord(scissor->w);
            scissor_y2 = scissor_y1 - get_coord(scissor->h);
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(LINE, line)
            float x1 = get_coord(line->begin.x);
            float y1 = 1.f - get_coord(line->begin.y);
            float x2 = get_coord(line->end.x);
            float y2 = 1.f - get_coord(line->end.y);
            struct Color color = get_color(line->color);
            float thickness = get_coord(line->line_thickness);
            render_line_screen(x1, y1, x2, y2, thickness, color, context);
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(CURVE, curve)
            printf("Render curve \n");
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(RECT, rect)
            struct Color color = get_color(rect->color);
            float thickness = get_coord(rect->line_thickness);
            float x1 = get_coord(rect->x)-0.5f*thickness;
            float y1 = 1.f - get_coord(rect->y) + 0.5f*thickness;
            float w = get_coord(rect->w) + thickness;
            float h = -get_coord(rect->h) - thickness;
            struct Matrix3 m = {
                w, 0.f, 0.f,
                0.f,   h, 0.f,
                x1,  y1, 1.f,
            };
            float r = get_coord(rect->rounding);
            if(r < 1e-3f) r = 1e-3f;
            struct Vec2 radius = {r/w, -r/h};
            int outline = 1;
            struct Vec2 thickness2 = {thickness/w, -thickness/h};
            struct ShaderUniform uniforms[] = {
                {"color", SHADER_UNIFORM_VEC3, 1, &color.r},
                {"radius", SHADER_UNIFORM_VEC2, 1, &radius},
                {"outline", SHADER_UNIFORM_INT, 1, &outline},
                {"thickness", SHADER_UNIFORM_VEC2, 1, &thickness2},
            };
            int num_uniforms = sizeof(uniforms)/sizeof(*uniforms);
            render_quad(gui->rect_shader, m, uniforms, num_uniforms, context);
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(RECT_FILLED, rect_filled)
            struct Color color = get_color(rect_filled->color);
            float x1 = get_coord(rect_filled->x);
            float y1 = 1.f - get_coord(rect_filled->y);
            float x2 = x1 + get_coord(rect_filled->w);
            float y2 = y1 - get_coord(rect_filled->h);
            //SCISSOR(x1, x) SCISSOR(x2, x) SCISSOR(y1, y) SCISSOR(y2, y)
            float w = x2 - x1;
            float h = y2 - y1;
            struct Matrix3 m = {
                  w, 0.f, 0.f,
                0.f,   h, 0.f,
                 x1,  y1, 1.f,
            };
            float r = get_coord(rect_filled->rounding);
            struct Vec2 radius = {r/w, -r/h};
            int outline = 0;
            struct ShaderUniform uniforms[] = {
                {"color", SHADER_UNIFORM_VEC3, 1, &color.r},
                {"radius", SHADER_UNIFORM_VEC2, 1, &radius},
                {"outline", SHADER_UNIFORM_INT, 1, &outline},
            };
            int num_uniforms = sizeof(uniforms)/sizeof(*uniforms);
            render_quad(gui->rect_shader, m, uniforms, num_uniforms, context);
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(RECT_MULTI_COLOR, rect_multi_color)
            struct Color colors[4];
            colors[0] = get_color(rect_multi_color->right);
            colors[1] = get_color(rect_multi_color->bottom);
            colors[2] = get_color(rect_multi_color->left);
            colors[3] = get_color(rect_multi_color->top);
            float x1 = get_coord(rect_multi_color->x);
            float y1 = 1.f - get_coord(rect_multi_color->y);
            float x2 = x1 + get_coord(rect_multi_color->w);
            float y2 = y1 - get_coord(rect_multi_color->h);
            float w = x2 - x1;
            float h = y2 - y1;
            struct Matrix3 m = {
                w, 0.f, 0.f,
                0.f,   h, 0.f,
                x1,  y1, 1.f,
            };
            struct ShaderUniform uniforms[] = {
                {"colors", SHADER_UNIFORM_VEC4, 4, colors},
            };
            int num_uniforms = sizeof(uniforms)/sizeof(*uniforms);
            render_quad(gui->color_rect_shader, m, uniforms, num_uniforms, context);
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(CIRCLE, circle)
            printf("Render circle\n");
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(CIRCLE_FILLED, circle_filled)
                struct Color color = get_color(circle_filled->color);
                float x1 = get_coord(circle_filled->x);
                float y1 = 1.f - get_coord(circle_filled->y);
                float x2 = x1 + get_coord(circle_filled->w);
                float y2 = y1 - get_coord(circle_filled->h);
                float w = x2 - x1;
                float h = y2 - y1;
                struct Matrix3 m = {
                    w, 0.f, 0.f,
                    0.f,   h, 0.f,
                    x1,  y1, 1.f,
                };
                struct ShaderUniform uniforms[] = {
                    {"color", SHADER_UNIFORM_VEC4, 1, &color.r},
                };
                int num_uniforms = sizeof(uniforms)/sizeof(*uniforms);
                render_quad(gui->circle_shader, m, uniforms, num_uniforms, context);
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(ARC, arc)
            printf("Render ??? \n");
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(ARC_FILLED, arc_filled)
            printf("Render ??? \n");
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(TRIANGLE, triangle)
            printf("Render ??? \n");
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(TRIANGLE_FILLED, triangle_filled)
            //printf("Render triangle filled \n");
            float x1 = get_coord(triangle_filled->a.x);
            float y1 = get_coord(triangle_filled->a.y);
            float x2 = get_coord(triangle_filled->b.x);
            float y2 = get_coord(triangle_filled->b.y);
            float x3 = get_coord(triangle_filled->c.x);
            float y3 = get_coord(triangle_filled->c.y);
            float x_min = x1 < x2 ? x1 : x2; x_min = x_min < x3 ? x_min : x3;
            float y_min = y1 < y2 ? y1 : y2; y_min = y_min < y3 ? y_min : y3;
            float x_max = x1 > x2 ? x1 : x2; x_max = x_max > x3 ? x_max : x3;
            float y_max = y1 > y2 ? y1 : y2; y_max = y_max > y3 ? y_max : y3;
            float x = x_min;
            float y = y_min;
            float w = x_max - x_min;
            float h = y_max - y_min;
            struct Matrix3 m = {
                w, 0.f, 0.f,
                0.f,   -h, 0.f,
                x,  1.f-y, 1.f,
            };
            struct Vec2 a = {(x1 - x_min)/w, (y1 - y_min)/h};
            struct Vec2 b = {(x2 - x_min)/w, (y2 - y_min)/h};
            struct Vec2 c = {(x3 - x_min)/w, (y3 - y_min)/h};
            struct Color color = get_color(triangle_filled->color);
            struct ShaderUniform uniforms[] = {
                {"color", SHADER_UNIFORM_VEC3, 1, &color.r},
                {"a",     SHADER_UNIFORM_VEC2, 1, &a},
                {"b",     SHADER_UNIFORM_VEC2, 1, &b},
                {"c",     SHADER_UNIFORM_VEC2, 1, &c},
            };
            int num_uniforms = sizeof(uniforms)/sizeof(*uniforms);
            render_quad(gui->triangle_shader, m, uniforms, num_uniforms, context);
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(POLYGON, polygon)
            printf("Render ??? \n");
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(POLYGON_FILLED, polygon_filled)
            printf("Render ??? \n");
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(POLYLINE, polyline)
            printf("Render polyline \n");
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(TEXT, text)
                struct GuiFont *font = text->font->userdata.ptr;
                float x = get_coord(text->x);
                float y = 1.f-get_coord(text->y);
                float h = get_font_height(font->font, context->data);
                struct Color color = get_color(text->foreground);
                //printf("rendering string \"%s\"\n", text->string);
                //render_line_screen(x, y, x+0.01f, y, 0.01f, color, context);
                y -= h/(float)reference_resolution;
                render_string_screen_n(text->string, text->length, font->font, &x,
                    &y, color, context);
        NK_COMMAND_CASE_END
         NK_COMMAND_CASE_BEGIN(IMAGE, image)
            printf("Render ??? \n");
        NK_COMMAND_CASE_END
        NK_COMMAND_CASE_BEGIN(CUSTOM, custom)
            printf("Render ??? \n");
        NK_COMMAND_CASE_END
        default:
            break;
        }
    }
    nk_clear(&gui->ctx);
}

void gui_fps_section(struct GuiContext *gui, struct GameData *data)
{
    struct nk_context *ctx = &gui->ctx;
    char buffer[32];
    int num_ticks = 64;
    int fps = 0;
    int mean_ticks = get_mean_ticks(num_ticks,data);
    if(mean_ticks > 0){
        fps = (ticks_per_second + mean_ticks/2)/mean_ticks;
    }
    sprintf(buffer, "fps: %d", fps);
    if(nk_tree_push(ctx, NK_TREE_TAB, buffer, NK_MINIMIZED)){
        assert(get_tick_length(data) >= num_ticks);
        nk_layout_row_dynamic(ctx, 150, 1);
        nk_chart_begin(ctx, NK_CHART_COLUMN, num_ticks, 0.f, 32.f);
        nk_chart_add_slot(ctx, NK_CHART_LINES, 1, 0.f, 32.f);
        nk_chart_push_slot(ctx, 1000.f/60.f, 1);
        nk_chart_push_slot(ctx, 1000.f/60.f, 1);
        int hovering_index = -1;
        for (int i = 0; i < num_ticks; ++i)
        {
            int ticks = get_ticks(i,data);
            float value = (float)ticks*inv_ticks_per_second*1000.f;
            nk_flags res = nk_chart_push(ctx, value);
            if(res & NK_CHART_HOVERING){
                hovering_index = i;
            }
        }
        nk_chart_end(ctx);
        /*
        if(hovering_index != -1){
            int ticks = get_ticks(hovering_index,data);
            float value = (float)ticks*inv_ticks_per_second*1000.f;
            nk_tooltipf(ctx, "blaj", value);
        }
         */
        nk_tree_pop(ctx);
    }
}

struct Color gui_color_picker(struct Color col, struct GuiContext *gui,
    struct GameData *data)
{
    struct nk_context *ctx = &gui->ctx;
    struct nk_colorf combo_col = {col.r, col.g, col.b, col.a};
    if (nk_combo_begin_color(ctx, nk_rgb_cf(combo_col), nk_vec2(200,400))) {
        enum color_mode {COL_RGB, COL_HSV};
        static int col_mode = COL_RGB;
        nk_layout_row_dynamic(ctx, 120, 1);
        combo_col = nk_color_picker(ctx, combo_col, NK_RGBA);
        
        nk_layout_row_dynamic(ctx, 25, 2);
        col_mode = nk_option_label(ctx, "RGB", col_mode == COL_RGB) ? COL_RGB : col_mode;
        col_mode = nk_option_label(ctx, "HSV", col_mode == COL_HSV) ? COL_HSV : col_mode;
        
        nk_layout_row_dynamic(ctx, 25, 1);
        if (col_mode == COL_RGB) {
            combo_col.r = nk_propertyf(ctx, "#R:", 0, combo_col.r, 1.0f, 0.01f,0.005f);
            combo_col.g = nk_propertyf(ctx, "#G:", 0, combo_col.g, 1.0f, 0.01f,0.005f);
            combo_col.b = nk_propertyf(ctx, "#B:", 0, combo_col.b, 1.0f, 0.01f,0.005f);
            combo_col.a = nk_propertyf(ctx, "#A:", 0, combo_col.a, 1.0f, 0.01f,0.005f);
        } else {
            float hsva[4];
            nk_colorf_hsva_fv(hsva, combo_col);
            hsva[0] = nk_propertyf(ctx, "#H:", 0, hsva[0], 1.0f, 0.01f,0.05f);
            hsva[1] = nk_propertyf(ctx, "#S:", 0, hsva[1], 1.0f, 0.01f,0.05f);
            hsva[2] = nk_propertyf(ctx, "#V:", 0, hsva[2], 1.0f, 0.01f,0.05f);
            hsva[3] = nk_propertyf(ctx, "#A:", 0, hsva[3], 1.0f, 0.01f,0.05f);
            combo_col = nk_hsva_colorfv(hsva);
        }
        nk_combo_end(ctx);
    }
    struct Color ret = {combo_col.r, combo_col.g, combo_col.b, combo_col.a};
    return ret;
}
