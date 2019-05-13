from matplotlib import pyplot as plt
import numpy as np
from utils import log
logger = log.get_console_logger()


def get_renderer(fig):
    try:
        return fig.canvas.get_renderer()
    except AttributeError:
        return fig.canvas.renderer


def get_midpoint(bbox):
    cx = (bbox.x0+bbox.x1)/2
    cy = (bbox.y0+bbox.y1)/2
    return cx, cy


def overlap_bbox_and_point(bbox, xp, yp):
    """Given a bbox that contains a given point, return the (x, y) displacement
    necessary to make the bbox not overlap the point."""
    cx, cy = get_midpoint(bbox)

    dir_x = np.sign(cx-xp)
    dir_y = np.sign(cy-yp)

    if dir_x == -1:
        dx = xp - bbox.xmax
    elif dir_x == 1:
        dx = xp - bbox.xmin
    else:
        dx = 0

    if dir_y == -1:
        dy = yp - bbox.ymax
    elif dir_y == 1:
        dy = yp - bbox.ymin
    else:
        dy = 0
    return dx, dy


def get_text_position(text, ax=None):
    ax = ax or plt.gca()
    x, y = text.get_position()
    return (ax.xaxis.convert_units(x),
            ax.yaxis.convert_units(y))


def get_bboxes(objs, renderer=None, expand=(1.0, 1.0), ax=None):
    if ax is None:
        ax = plt.gca()
    if renderer is None:
        renderer = get_renderer(ax.get_figure())
    return [i.get_window_extent(renderer).expanded(*expand).transformed(ax.\
                                          transData.inverted()) for i in objs]


def get_points_inside_bbox(x, y, bbox):
    """Return the indices of points inside the given bbox."""
    x1, y1, x2, y2 = bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax
    x_in = np.logical_and(x>x1, x<x2)
    y_in = np.logical_and(y>y1, y<y2)
    return np.where(x_in & y_in)[0]


def move_texts(texts, delta_x, delta_y, bboxes=None, renderer=None, ax=None, constrain_to_axis_limits=True):
    if len(texts) != len(delta_x) or len(delta_x) != len(delta_y):
        raise ValueError("Input arguments texts, delta_x and delta_y must have the same size.")

    if ax is None:
        ax = plt.gca()
    if bboxes is None and constrain_to_axis_limits:
        if renderer is None:
            r = get_renderer(ax.get_figure())
        else:
            r = renderer
        bboxes = get_bboxes(texts, r, (1, 1), ax=ax)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    dx_arr = []
    dy_arr = []

    for i, (text, dx, dy) in enumerate(zip(texts, delta_x, delta_y)):
        if constrain_to_axis_limits:
            bbox = bboxes[i]
            x1, y1, x2, y2 = bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymax
            if x1 + dx < xmin:
                dx = 0
            if x2 + dx > xmax:
                dx = 0
            if y1 + dy < ymin:
                dy = 0
            if y2 + dy > ymax:
                dy = 0

        x, y = get_text_position(text, ax=ax)
        dx_arr.append(dx)
        dy_arr.append(dy)
        newx = x + dx
        newy = y + dy
        text.set_position((newx, newy))

    return np.array(dx_arr), np.array(dy_arr)


def rearrange_text_radially(texts, alpha=0.5, ax=None, renderer=None, move=True, origin=(0, 0), constrain_to_ax_limits=True):
    """

    :param texts:
    :param alpha: Weighting parameter balancing repulsion between text objects and the radial outward force.
    Must be in [0, 1], where 0 means 'no radial force' and 1 means 'all radial force'
    :param ax:
    :param renderer:
    :param move:
    :param origin:
    :return:
    """
    origin = np.array(origin)
    if ax is None:
        ax = plt.gca()
    if renderer is None:
        r = get_renderer(ax.get_figure())
    else:
        r = renderer
    bboxes = get_bboxes(texts, r, ax=ax)
    xmins = [bbox.xmin for bbox in bboxes]
    xmaxs = [bbox.xmax for bbox in bboxes]
    ymaxs = [bbox.ymax for bbox in bboxes]
    ymins = [bbox.ymin for bbox in bboxes]

    overlaps_x = np.zeros((len(bboxes), len(bboxes)))
    overlaps_y = np.zeros_like(overlaps_x)
    overlap_directions_x = np.zeros_like(overlaps_x)
    overlap_directions_y = np.zeros_like(overlaps_y)
    radial_x = np.zeros_like(overlaps_x)
    radial_y = np.zeros_like(overlaps_y)
    for i, bbox1 in enumerate(bboxes):
        overlaps = get_points_inside_bbox(xmins*2+xmaxs*2, (ymins+ymaxs)*2,
                                             bbox1) % len(bboxes)
        overlaps = np.unique(overlaps)
        for j in overlaps:
            bbox2 = bboxes[j]
            x, y = bbox1.intersection(bbox1, bbox2).size
            overlaps_x[i, j] = x
            overlaps_y[i, j] = y
            direction = np.sign(bbox1.extents - bbox2.extents)[:2]
            overlap_directions_x[i, j] = direction[0]
            overlap_directions_y[i, j] = direction[1]

            # normed radial vector
            rn = np.array(bbox1.extents[:2]) - origin
            rn /= (rn ** 2).sum() ** 0.5
            # project shift vector onto this, ensuring we're pointing outwards (take the absolute value)
            proj = np.abs(rn.dot([x, y])) * rn
            radial_x[i, j] = proj[0]
            radial_y[i, j] = proj[1]

    move_x = (1 - alpha) * overlaps_x*overlap_directions_x + alpha * radial_x
    move_y = (1 - alpha) * overlaps_y*overlap_directions_y + alpha * radial_y

    delta_x = move_x.sum(axis=1)
    delta_y = move_y.sum(axis=1)

    q = np.sum(overlaps_x), np.sum(overlaps_y)
    if move:
        delta_x, delta_y = move_texts(texts, delta_x, delta_y, bboxes, ax=ax, constrain_to_axis_limits=constrain_to_ax_limits)

    return delta_x, delta_y, q


def adjust_text_radial_plus_repulsion(
    texts,
    alpha=0.5,
    ax=None,
    add_line=True,
    origin=(0, 0),
    n_iter_max=15,
    tol=1e-4,
    lineprops=None,
    only_draw_if_non_intersecting=True,
    constrain_to_ax_limits=False
):
    # TODO: would like to add repulsion between text and points, too
    """

    :param texts:
    :param alpha:
    :param add_line:
    :param origin:
    :param n_iter_max:
    :param tol: The value of the L1 norm of all overlap vectors at which we consider the job done. This will depend upon
    the scale of the plot. TODO: could select it automatically based on a heuristic?
    :return:
    """
    if add_line and lineprops is None:
        lineprops = {
            'color': 'black',
            'linewidth': 1.,
        }

    if ax is None:
        ax = plt.gca()

    if add_line and only_draw_if_non_intersecting:
        r = get_renderer(ax.get_figure())
        bboxes_orig = get_bboxes(texts, r, ax=ax)

    # record original positions
    orig_pos = [t.get_position() for t in texts]
    dx = None
    dy = None
    converged = False

    for i in range(n_iter_max):
        this_dx, this_dy, q = rearrange_text_radially(texts, alpha=alpha, origin=origin, constrain_to_ax_limits=constrain_to_ax_limits)

        if dx is None:
            dx = this_dx
        else:
            dx += this_dx
        if dy is None:
            dy = this_dy
        else:
            dy += this_dy

        if sum(q) < tol:
            logger.info("Rearrange text converged after %d iterations.", (i + 1))
            converged = True
            break

    if not converged:
        logger.warning("Failed to converge after %d iterations", n_iter_max)

    if add_line:
        if only_draw_if_non_intersecting:
            r = get_renderer(ax.get_figure())
            bboxes = get_bboxes(texts, r, ax=ax)

        for i, (xy, ddx, ddy) in enumerate(zip(orig_pos, dx, dy)):
            # check (1): is the line long enough to be worth drawing?
            if np.abs(np.array([ddx, ddy])).sum() < tol:
                continue

            new_x = xy[0] + ddx
            new_y = xy[1] + ddy

            # check (2): is the new location still intersecting the original position?
            if only_draw_if_non_intersecting:
                this_bbox = bboxes[i]
                this_bbox_orig = bboxes_orig[i]
                if this_bbox.overlaps(this_bbox_orig):
                    continue
                ## FIXME: this isn't checking intersection
                # if this_bbox.contains(*xy):
                #     continue

            ax.plot(
                [xy[0], new_x],
                [xy[1], new_y],
                **lineprops
            )

    return dx, dy, (i + 1)





