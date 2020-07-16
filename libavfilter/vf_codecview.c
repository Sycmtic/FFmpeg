/*
 * Copyright (c) 2002-2004 Michael Niedermayer <michaelni@gmx.at>
 * Copyright (c) 2014 Clément Bœsch <u pkh me>
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/**
 * @file
 * Codec debug viewer filter.
 *
 * All the MV drawing code from Michael Niedermayer is extracted from
 * libavcodec/mpegvideo.c.
 *
 * TODO: segmentation
 */

#include "libavutil/imgutils.h"
#include "libavutil/motion_vector.h"
#include "libavutil/opt.h"
#include "avfilter.h"
#include "internal.h"
#include "libavutil/video_enc_params.h"

#define MV_P_FOR  (1<<0)
#define MV_B_FOR  (1<<1)
#define MV_B_BACK (1<<2)
#define MV_TYPE_FOR  (1<<0)
#define MV_TYPE_BACK (1<<1)
#define FRAME_TYPE_I (1<<0)
#define FRAME_TYPE_P (1<<1)
#define FRAME_TYPE_B (1<<2)

typedef struct CodecViewContext {
    const AVClass *class;
    unsigned mv;
    unsigned frame_type;
    unsigned mv_type;
    int hsub, vsub;
    int qp;
    int chroma_qp;
    int dc_qp;
    int bs;
} CodecViewContext;

#define OFFSET(x) offsetof(CodecViewContext, x)
#define FLAGS AV_OPT_FLAG_FILTERING_PARAM|AV_OPT_FLAG_VIDEO_PARAM
#define CONST(name, help, val, unit) { name, help, 0, AV_OPT_TYPE_CONST, {.i64=val}, 0, 0, FLAGS, unit }

static const AVOption codecview_options[] = {
    { "mv", "set motion vectors to visualize", OFFSET(mv), AV_OPT_TYPE_FLAGS, {.i64=0}, 0, INT_MAX, FLAGS, "mv" },
        CONST("pf", "forward predicted MVs of P-frames",  MV_P_FOR,  "mv"),
        CONST("bf", "forward predicted MVs of B-frames",  MV_B_FOR,  "mv"),
        CONST("bb", "backward predicted MVs of B-frames", MV_B_BACK, "mv"),
    { "qp", NULL, OFFSET(qp), AV_OPT_TYPE_BOOL, {.i64=0}, 0, 1, .flags = FLAGS },
    { "chroma_qp", NULL, OFFSET(chroma_qp), AV_OPT_TYPE_BOOL, {.i64=0}, 0, 1, .flags = FLAGS },
    { "dc_qp", NULL, OFFSET(dc_qp), AV_OPT_TYPE_BOOL, {.i64=0}, 0, 1, .flags = FLAGS },
    { "bs", "set block structure to visualize", OFFSET(bs), AV_OPT_TYPE_BOOL, {.i64=0}, 0, 1, .flags = FLAGS },
    { "mv_type", "set motion vectors type", OFFSET(mv_type), AV_OPT_TYPE_FLAGS, {.i64=0}, 0, INT_MAX, FLAGS, "mv_type" },
    { "mvt",     "set motion vectors type", OFFSET(mv_type), AV_OPT_TYPE_FLAGS, {.i64=0}, 0, INT_MAX, FLAGS, "mv_type" },
        CONST("fp", "forward predicted MVs",  MV_TYPE_FOR,  "mv_type"),
        CONST("bp", "backward predicted MVs", MV_TYPE_BACK, "mv_type"),
    { "frame_type", "set frame types to visualize motion vectors of", OFFSET(frame_type), AV_OPT_TYPE_FLAGS, {.i64=0}, 0, INT_MAX, FLAGS, "frame_type" },
    { "ft",         "set frame types to visualize motion vectors of", OFFSET(frame_type), AV_OPT_TYPE_FLAGS, {.i64=0}, 0, INT_MAX, FLAGS, "frame_type" },
        CONST("if", "I-frames", FRAME_TYPE_I, "frame_type"),
        CONST("pf", "P-frames", FRAME_TYPE_P, "frame_type"),
        CONST("bf", "B-frames", FRAME_TYPE_B, "frame_type"),
    { NULL }
};

AVFILTER_DEFINE_CLASS(codecview);

static int query_formats(AVFilterContext *ctx)
{
    // TODO: we can probably add way more pixel formats without any other
    // changes; anything with 8-bit luma in first plane should be working
    static const enum AVPixelFormat pix_fmts[] = {AV_PIX_FMT_YUV420P, AV_PIX_FMT_NONE};
    AVFilterFormats *fmts_list = ff_make_format_list(pix_fmts);
    if (!fmts_list)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, fmts_list);
}

static int clip_line(int *sx, int *sy, int *ex, int *ey, int maxx)
{
    if(*sx > *ex)
        return clip_line(ex, ey, sx, sy, maxx);

    if (*sx < 0) {
        if (*ex < 0)
            return 1;
        *sy = *ey + (*sy - *ey) * (int64_t)*ex / (*ex - *sx);
        *sx = 0;
    }

    if (*ex > maxx) {
        if (*sx > maxx)
            return 1;
        *ey = *sy + (*ey - *sy) * (int64_t)(maxx - *sx) / (*ex - *sx);
        *ex = maxx;
    }
    return 0;
}

/**
 * Draw a line from (ex, ey) -> (sx, sy).
 * @param w width of the image
 * @param h height of the image
 * @param stride stride/linesize of the image
 * @param color color of the arrow
 */
static void draw_line(uint8_t *buf, int sx, int sy, int ex, int ey,
                      int w, int h, int stride, int color)
{
    int x, y, fr, f;

    if (clip_line(&sx, &sy, &ex, &ey, w - 1))
        return;
    if (clip_line(&sy, &sx, &ey, &ex, h - 1))
        return;

    sx = av_clip(sx, 0, w - 1);
    sy = av_clip(sy, 0, h - 1);
    ex = av_clip(ex, 0, w - 1);
    ey = av_clip(ey, 0, h - 1);

    buf[sy * stride + sx] += color;

    if (FFABS(ex - sx) > FFABS(ey - sy)) {
        if (sx > ex) {
            FFSWAP(int, sx, ex);
            FFSWAP(int, sy, ey);
        }
        buf += sx + sy * stride;
        ex  -= sx;
        f    = ((ey - sy) << 16) / ex;
        for (x = 0; x <= ex; x++) {
            y  = (x * f) >> 16;
            fr = (x * f) & 0xFFFF;
                   buf[ y      * stride + x] += (color * (0x10000 - fr)) >> 16;
            if(fr) buf[(y + 1) * stride + x] += (color *            fr ) >> 16;
        }
    } else {
        if (sy > ey) {
            FFSWAP(int, sx, ex);
            FFSWAP(int, sy, ey);
        }
        buf += sx + sy * stride;
        ey  -= sy;
        if (ey)
            f = ((ex - sx) << 16) / ey;
        else
            f = 0;
        for(y= 0; y <= ey; y++){
            x  = (y*f) >> 16;
            fr = (y*f) & 0xFFFF;
                   buf[y * stride + x    ] += (color * (0x10000 - fr)) >> 16;
            if(fr) buf[y * stride + x + 1] += (color *            fr ) >> 16;
        }
    }
}

/**
 * Draw an arrow from (ex, ey) -> (sx, sy).
 * @param w width of the image
 * @param h height of the image
 * @param stride stride/linesize of the image
 * @param color color of the arrow
 */
static void draw_arrow(uint8_t *buf, int sx, int sy, int ex,
                       int ey, int w, int h, int stride, int color, int tail, int direction)
{
    int dx,dy;

    if (direction) {
        FFSWAP(int, sx, ex);
        FFSWAP(int, sy, ey);
    }

    sx = av_clip(sx, -100, w + 100);
    sy = av_clip(sy, -100, h + 100);
    ex = av_clip(ex, -100, w + 100);
    ey = av_clip(ey, -100, h + 100);

    dx = ex - sx;
    dy = ey - sy;

    if (dx * dx + dy * dy > 3 * 3) {
        int rx =  dx + dy;
        int ry = -dx + dy;
        int length = sqrt((rx * rx + ry * ry) << 8);

        // FIXME subpixel accuracy
        rx = ROUNDED_DIV(rx * 3 << 4, length);
        ry = ROUNDED_DIV(ry * 3 << 4, length);

        if (tail) {
            rx = -rx;
            ry = -ry;
        }

        draw_line(buf, sx, sy, sx + rx, sy + ry, w, h, stride, color);
        draw_line(buf, sx, sy, sx - ry, sy + rx, w, h, stride, color);
    }
    draw_line(buf, sx, sy, ex, ey, w, h, stride, color);
}

static int qp_color_calculate(int qp, enum AVVideoEncParamsType type) {
    return type == AV_VIDEO_ENC_PARAMS_H264 ? qp * 128 / 31 : qp;
}

static void get_block_color(AVVideoEncParams *par, AVVideoBlockParams *b, CodecViewContext *s, enum AVColorRange color_range, int *cu, int *cv)
{
    const int plane_qp_cu_index = s->chroma_qp ? 1 : 0;
    const int plane_qp_cv_index = s->chroma_qp ? 2 : 0;
    const int ac_dc_index = s->dc_qp ? 0 : 1;
    *cu = qp_color_calculate(par->qp + par->delta_qp[plane_qp_cu_index][ac_dc_index] + b->delta_qp, par->type);
    *cv = qp_color_calculate(par->qp + par->delta_qp[plane_qp_cv_index][ac_dc_index] + b->delta_qp, par->type);
    if (color_range == AVCOL_RANGE_MPEG) {
        // map jpeg color range(0-255) to mpeg color range(16-235)
        *cu = av_rescale(*cu, 73, 85) + 16;
        *cv = av_rescale(*cv, 73, 85) + 16;
    }
}

static void color_block(AVFrame *frame, CodecViewContext *s, const int src_x, const int src_y, const int b_w, const int b_h, const int cu, const int cv) {
    const int w = AV_CEIL_RSHIFT(frame->width, s->hsub);
    const int h = AV_CEIL_RSHIFT(frame->height, s->vsub);
    const int lzu = frame->linesize[1];
    const int lzv = frame->linesize[2];

    const int plane_src_x = src_x >> s->hsub;
    const int plane_src_y = src_y >> s->vsub;
    const int plane_b_w = b_w >> s->hsub;
    const int plane_b_h = b_h >> s->vsub;
    uint8_t *pu = frame->data[1] + plane_src_y * lzu;
    uint8_t *pv = frame->data[2] + plane_src_y * lzv;

    for (int y = plane_src_y; y < plane_src_y + plane_b_h; y++) {
        for (int x = plane_src_x; x < plane_src_x + plane_b_w; x++) {
            if (x >= w)
                break;
            pu[x] = cu;
            pv[x] = cv;
        }
        if (y >= h)
            break;
        pu += lzu;
        pv += lzv;
    }
}

static void draw_block_border(AVFrame *frame, AVVideoBlockParams *b)
{
    const int lzy = frame->linesize[0];
    uint8_t *py = frame->data[0] + b->src_y * lzy;

    for (int x = b->src_x; x < b->src_x + b->w; x++) {
        if (x >= frame->width)
            break;
        py[x] = py[x] * 3 / 4;
    }
    for (int y = b->src_y; y < b->src_y + b->h; y++) {
        if (y >= frame->height)
            break;
        py[b->src_x] = py[b->src_x] * 3 / 4;
        py[b->src_x + b->w - 1] = py[b->src_x + b->w - 1] * 3 / 4;
        py += lzy;
    }
    for (int x = b->src_x; x < b->src_x + b->w; x++) {
        if (x >= frame->width)
            break;
        py[x] = py[x] * 3 / 4;
    }
}

static int filter_frame(AVFilterLink *inlink, AVFrame *frame)
{
    AVFilterContext *ctx = inlink->dst;
    CodecViewContext *s = ctx->priv;
    AVFilterLink *outlink = ctx->outputs[0];

    AVFrameSideData *sd = av_frame_get_side_data(frame, AV_FRAME_DATA_VIDEO_ENC_PARAMS);

    if (s->qp) {
        int qstride, qp_type;
        int8_t *qp_table = av_frame_get_qp_table(frame, &qstride, &qp_type);

        if (qp_table) {
            int x, y;
            const int w = AV_CEIL_RSHIFT(frame->width, s->hsub);
            const int h = AV_CEIL_RSHIFT(frame->height, s->vsub);
            uint8_t *pu = frame->data[1];
            uint8_t *pv = frame->data[2];
            const int lzu = frame->linesize[1];
            const int lzv = frame->linesize[2];

            for (y = 0; y < h; y++) {
                for (x = 0; x < w; x++) {
                    const int qp = ff_norm_qscale(qp_table[(y >> 3) * qstride + (x >> 3)], qp_type) * 128 / 31;
                    pu[x] = pv[x] = qp;
                }
                pu += lzu;
                pv += lzv;
            }
        }

        if (sd) {
            AVVideoEncParams *par = (AVVideoEncParams*)sd->data;

            if (par->nb_blocks) {
                for (int i = 0; i < par->nb_blocks; i++) {
                    AVVideoBlockParams *b = av_video_enc_params_block(par, i);
                    int cu, cv;
                    get_block_color(par, b, s, frame->color_range, &cu, &cv);
                    color_block(frame, s, b->src_x, b->src_y, b->w, b->h, cu, cv);
                }
            } else {
                const c = qp_color_calculate(par->qp, par->type);
                color_block(frame, s, 0, 0, frame->width, frame->height, c, c);
            }
        }
    }

    if (s->bs) {
        if (sd) {
            AVVideoEncParams *par = (AVVideoEncParams*)sd->data;

            if (par->nb_blocks) {
                for (int i = 0; i < par->nb_blocks; i++) {
                    AVVideoBlockParams *b = av_video_enc_params_block(par, i);
                    draw_block_border(frame, b);
                }
            }
        }
    }

    if (s->mv || s->mv_type) {
        AVFrameSideData *sd = av_frame_get_side_data(frame, AV_FRAME_DATA_MOTION_VECTORS);
        if (sd) {
            int i;
            const AVMotionVector *mvs = (const AVMotionVector *)sd->data;
            const int is_iframe = (s->frame_type & FRAME_TYPE_I) && frame->pict_type == AV_PICTURE_TYPE_I;
            const int is_pframe = (s->frame_type & FRAME_TYPE_P) && frame->pict_type == AV_PICTURE_TYPE_P;
            const int is_bframe = (s->frame_type & FRAME_TYPE_B) && frame->pict_type == AV_PICTURE_TYPE_B;

            for (i = 0; i < sd->size / sizeof(*mvs); i++) {
                const AVMotionVector *mv = &mvs[i];
                const int direction = mv->source > 0;

                if (s->mv_type) {
                    const int is_fp = direction == 0 && (s->mv_type & MV_TYPE_FOR);
                    const int is_bp = direction == 1 && (s->mv_type & MV_TYPE_BACK);

                    if ((!s->frame_type && (is_fp || is_bp)) ||
                        is_iframe && is_fp || is_iframe && is_bp ||
                        is_pframe && is_fp ||
                        is_bframe && is_fp || is_bframe && is_bp)
                        draw_arrow(frame->data[0], mv->dst_x, mv->dst_y, mv->src_x, mv->src_y,
                                   frame->width, frame->height, frame->linesize[0],
                                   100, 0, direction);
                } else if (s->mv)
                    if ((direction == 0 && (s->mv & MV_P_FOR)  && frame->pict_type == AV_PICTURE_TYPE_P) ||
                        (direction == 0 && (s->mv & MV_B_FOR)  && frame->pict_type == AV_PICTURE_TYPE_B) ||
                        (direction == 1 && (s->mv & MV_B_BACK) && frame->pict_type == AV_PICTURE_TYPE_B))
                        draw_arrow(frame->data[0], mv->dst_x, mv->dst_y, mv->src_x, mv->src_y,
                                   frame->width, frame->height, frame->linesize[0],
                                   100, 0, direction);
            }
        }
    }

    return ff_filter_frame(outlink, frame);
}

static int config_input(AVFilterLink *inlink)
{
    AVFilterContext *ctx = inlink->dst;
    CodecViewContext *s = ctx->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);

    s->hsub = desc->log2_chroma_w;
    s->vsub = desc->log2_chroma_h;
    return 0;
}

static const AVFilterPad codecview_inputs[] = {
    {
        .name           = "default",
        .type           = AVMEDIA_TYPE_VIDEO,
        .filter_frame   = filter_frame,
        .config_props   = config_input,
        .needs_writable = 1,
    },
    { NULL }
};

static const AVFilterPad codecview_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_codecview = {
    .name          = "codecview",
    .description   = NULL_IF_CONFIG_SMALL("Visualize information about some codecs."),
    .priv_size     = sizeof(CodecViewContext),
    .query_formats = query_formats,
    .inputs        = codecview_inputs,
    .outputs       = codecview_outputs,
    .priv_class    = &codecview_class,
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_GENERIC,
};
