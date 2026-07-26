import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from scipy.interpolate import RegularGridInterpolator


class InteractiveLineout:
    """
    Interactively extract a finite-width lineout from 2D gridded data.

    Interaction
    -----------
    Left-click and drag:
        Define the lineout centerline.

    Mouse wheel:
        Change the averaging width.

    Up / Down arrows:
        Increase / decrease the averaging width.

    Enter:
        Print the current lineout coordinates.

    Escape:
        Clear the current selection.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Existing axes containing the imshow, pcolormesh, contour plot, etc.

    x, y : 1D array-like
        Coordinates of the data.

        These may be either:
        - cell-center coordinates, with lengths nx and ny, or
        - cell-edge coordinates, with lengths nx + 1 and ny + 1.

    data : 2D array-like
        Data with shape (ny, nx), meaning data[j, i] corresponds to
        coordinate (x[i], y[j]).

    width : float, optional
        Initial physical averaging width in the same units as x and y.
        If omitted, approximately five local grid spacings are used.

    n_along : int or None, optional
        Number of samples along the selected centerline. If None, it is
        chosen automatically from the grid spacing.

    n_across : int or None, optional
        Number of samples across the averaging width. If None, it is
        chosen automatically from the grid spacing.

    interpolation : {"linear", "nearest"}, optional
        Interpolation method passed to RegularGridInterpolator.

    reduction : callable, optional
        Function used to combine samples across the strip.
        Typical choices are np.nanmean, np.nanmedian, or np.nanmax.

    profile_ax : matplotlib.axes.Axes or None, optional
        Existing axes on which to plot the lineout. If omitted, a new
        figure and axes are created.

    lineout_label : str, optional
        Vertical-axis label for the lineout plot.
    """

    def __init__(
        self,
        ax,
        x,
        y,
        data,
        width=None,
        n_along=None,
        n_across=None,
        interpolation="linear",
        reduction=np.nanmean,
        profile_ax=None,
        lineout_label="Value",
    ):
        self.ax = ax
        self.fig = ax.figure

        self.data = np.asarray(data, dtype=float)

        if self.data.ndim != 2:
            raise ValueError(
                f"data must be two-dimensional; got shape {self.data.shape}."
            )

        self.y, self.x = self._prepare_coordinates(
            x=np.asarray(x, dtype=float),
            y=np.asarray(y, dtype=float),
            data_shape=self.data.shape,
        )

        self.data, self.x, self.y = self._make_coordinates_ascending(
            self.data,
            self.x,
            self.y,
        )

        self.n_along_requested = n_along
        self.n_across_requested = n_across
        self.reduction = reduction

        self.dx = self._representative_spacing(self.x)
        self.dy = self._representative_spacing(self.y)

        self.max_along_samples = 2000
        self.max_across_samples = 200

        if width is None:
            width = 5.0 * min(self.dx,self.dy,)

        if width <= 0:
            raise ValueError("width must be positive.")

        self.width = float(width)

        self.interpolator = RegularGridInterpolator(
            (self.y, self.x),
            self.data,
            method=interpolation,
            bounds_error=False,
            fill_value=np.nan,
        )

        # Current selection and output
        self.start = None
        self.end = None
        self.distance = None
        self.profile = None
        self.x_line = None
        self.y_line = None
        self.sample_values = None

        # Mouse-drag state
        self._dragging = False

        # Artists drawn on the source axes
        self.centerline_artist, = self.ax.plot(
            [],
            [],
            linestyle="--",
            linewidth=1.5,
            marker="o",
            markersize=4,
            visible=False,
        )

        self.strip_artist = Polygon(
            np.zeros((4, 2)),
            closed=True,
            fill=False,
            linewidth=1.5,
            visible=False,
        )
        self.ax.add_patch(self.strip_artist)

        # Profile output figure
        if profile_ax is None:
            self.profile_fig, self.profile_ax = plt.subplots(
                figsize=(7, 4.5)
            )
        else:
            self.profile_ax = profile_ax
            self.profile_fig = profile_ax.figure

        self.profile_artist, = self.profile_ax.plot([], [])
        self.profile_ax.set_xlabel("Distance along line")
        self.profile_ax.set_ylabel(lineout_label)
        self.profile_ax.set_title("Averaged lineout")
        self.profile_ax.grid(alpha=0.3)

        # Status text
        self.status_artist = self.ax.text(
            0.02,
            0.98,
            self._status_message(),
            transform=self.ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            bbox={
                "facecolor": "white",
                "alpha": 0.75,
                "edgecolor": "none",
            },
        )

        # Keep callback IDs so the tool can later be disconnected cleanly.
        self._callback_ids = [
            self.fig.canvas.mpl_connect(
                "button_press_event",
                self._on_press,
            ),
            self.fig.canvas.mpl_connect(
                "motion_notify_event",
                self._on_motion,
            ),
            self.fig.canvas.mpl_connect(
                "button_release_event",
                self._on_release,
            ),
            self.fig.canvas.mpl_connect(
                "scroll_event",
                self._on_scroll,
            ),
            self.fig.canvas.mpl_connect(
                "key_press_event",
                self._on_key,
            ),
        ]

        self.fig.canvas.draw_idle()
        self.profile_fig.canvas.draw_idle()

    @staticmethod
    def _cell_centers(edges):
        """Convert cell-edge coordinates to cell-center coordinates."""
        return 0.5 * (edges[:-1] + edges[1:])

    @classmethod
    def _prepare_coordinates(cls, x, y, data_shape):
        """
        Accept center arrays or pcolormesh-style edge arrays.
        """
        ny, nx = data_shape

        if x.ndim != 1 or y.ndim != 1:
            raise ValueError(
                "x and y must be one-dimensional coordinate arrays."
            )

        if len(x) == nx + 1:
            x = cls._cell_centers(x)
        elif len(x) != nx:
            raise ValueError(
                f"x has length {len(x)}, but data has {nx} columns. "
                f"Expected {nx} centers or {nx + 1} edges."
            )

        if len(y) == ny + 1:
            y = cls._cell_centers(y)
        elif len(y) != ny:
            raise ValueError(
                f"y has length {len(y)}, but data has {ny} rows. "
                f"Expected {ny} centers or {ny + 1} edges."
            )

        return y, x

    @staticmethod
    def _make_coordinates_ascending(data, x, y):
        """
        RegularGridInterpolator accepts ascending or descending coordinates,
        but normalizing to ascending order simplifies the rest of the code.
        """
        if np.any(np.diff(x) == 0):
            raise ValueError("x coordinates must be strictly monotonic.")

        if np.any(np.diff(y) == 0):
            raise ValueError("y coordinates must be strictly monotonic.")

        if not (
            np.all(np.diff(x) > 0)
            or np.all(np.diff(x) < 0)
        ):
            raise ValueError("x coordinates must be monotonic.")

        if not (
            np.all(np.diff(y) > 0)
            or np.all(np.diff(y) < 0)
        ):
            raise ValueError("y coordinates must be monotonic.")

        if x[0] > x[-1]:
            x = x[::-1]
            data = data[:, ::-1]

        if y[0] > y[-1]:
            y = y[::-1]
            data = data[::-1, :]

        return data, x, y

    @staticmethod
    def _representative_spacing(coordinates):
        differences = np.abs(np.diff(coordinates))
        return float(np.median(differences))

    def _status_message(self):
        return (
            "Drag: select line\n"
            "Wheel or ↑/↓: change width\n"
            f"Width = {self.width:.4g}"
        )

    def _event_is_valid(self, event):
        return (
            event.inaxes is self.ax
            and event.xdata is not None
            and event.ydata is not None
        )

    def _on_press(self, event):
        if not self._event_is_valid(event):
            return

        if event.button != 1:
            return

        # Avoid conflicting with Matplotlib zoom and pan tools.
        toolbar = getattr(
            self.fig.canvas.manager,
            "toolbar",
            None,
        )

        if toolbar is not None and getattr(toolbar, "mode", ""):
            return

        self._dragging = True
        self.start = np.array(
            [event.xdata, event.ydata],
            dtype=float,
        )
        self.end = self.start.copy()

        self._update_selection()

    def _on_motion(self, event):
        if not self._dragging:
            return

        if not self._event_is_valid(event):
            return

        self.end = np.array(
            [event.xdata, event.ydata],
            dtype=float,
        )

        self._update_selection()

    def _on_release(self, event):
        if not self._dragging:
            return

        self._dragging = False

        if self._event_is_valid(event):
            self.end = np.array(
                [event.xdata, event.ydata],
                dtype=float,
            )

        self._update_selection()

    def _on_scroll(self, event):
        if event.inaxes is not self.ax:
            return

        if event.button == "up":
            self.width *= 1.15
        elif event.button == "down":
            self.width /= 1.15
        else:
            return

        self.width = max(
            self.width,
            0.1 * min(self.dx, self.dy),
        )

        self._update_selection()

    def _on_key(self, event):
        if event.key == "up":
            self.width *= 1.15
            self._update_selection()

        elif event.key == "down":
            self.width /= 1.15
            self.width = max(
                self.width,
                0.1 * min(self.dx, self.dy),
            )
            self._update_selection()

        elif event.key == "escape":
            self.clear()

        elif event.key == "enter":
            if self.distance is None:
                print("No lineout has been selected yet.")
            else:
                print(f"Start: {tuple(self.start)}")
                print(f"End:   {tuple(self.end)}")
                print(f"Width: {self.width}")
                print(f"Samples along line: {len(self.distance)}")

    def _selection_geometry(self):
        """
        Return line length, tangent unit vector, normal unit vector, and
        four rectangle corners.
        """
        delta = self.end - self.start
        length = float(np.hypot(delta[0], delta[1]))

        if length == 0:
            return None

        tangent = delta / length
        normal = np.array(
            [-tangent[1], tangent[0]]
        )

        half_width_vector = 0.5 * self.width * normal

        corners = np.array(
            [
                self.start + half_width_vector,
                self.end + half_width_vector,
                self.end - half_width_vector,
                self.start - half_width_vector,
            ]
        )

        return length, tangent, normal, corners

    def _choose_sampling_counts(self,length,tangent,normal,):

        """
        Choose safe sampling counts based on the grid resolution projected
        along and perpendicular to the selected line.
        """

        # Effective grid spacing along the selected line.
        spacing_along = np.sqrt(
            (tangent[0] * self.dx) ** 2
            + (tangent[1] * self.dy) ** 2
        )

        # Effective grid spacing perpendicular to the selected line.
        spacing_across = np.sqrt(
            (normal[0] * self.dx) ** 2
            + (normal[1] * self.dy) ** 2
        )

        spacing_along = max(
            spacing_along,
            np.finfo(float).eps,
        )

        spacing_across = max(
            spacing_across,
            np.finfo(float).eps,
        )

        if self.n_along_requested is None:
            n_along = (
                int(np.ceil(length / spacing_along))
                + 1
            )
        else:
            n_along = int(
                self.n_along_requested
            )

        if self.n_across_requested is None:
            n_across = (
                int(np.ceil(self.width / spacing_across))
                + 1
            )
        else:
            n_across = int(
                self.n_across_requested
            )

        # Keep the allocation safely bounded.
        n_along = np.clip(
            n_along,
            2,
            self.max_along_samples,
        )

        n_across = np.clip(
            n_across,
            1,
            self.max_across_samples,
        )

        return (
            int(n_along),
            int(n_across),
        )

    def _extract_profile(self, length, tangent, normal):
        n_along, n_across = self._choose_sampling_counts(length, tangent, normal)

        distance = np.linspace(
            0.0,
            length,
            n_along,
        )

        if n_across == 1:
            offsets = np.array([0.0])
        else:
            offsets = np.linspace(
                -0.5 * self.width,
                0.5 * self.width,
                n_across,
            )

        # Centerline points: shape (n_along, 2)
        center_points = (
            self.start[None, :]
            + distance[:, None] * tangent[None, :]
        )

        # Full strip coordinates: shape (n_along, n_across)
        x_samples = (
            center_points[:, 0, None]
            + offsets[None, :] * normal[0]
        )

        y_samples = (
            center_points[:, 1, None]
            + offsets[None, :] * normal[1]
        )

        # RegularGridInterpolator expects points ordered as (y, x).
        interpolation_points = np.column_stack(
            [
                y_samples.ravel(),
                x_samples.ravel(),
            ]
        )

        sampled = self.interpolator(
            interpolation_points
        ).reshape(n_along, n_across)

        with np.errstate(
            invalid="ignore",
            divide="ignore",
        ):
            profile = self.reduction(
                sampled,
                axis=1,
            )

        return (
            distance,
            profile,
            center_points[:, 0],
            center_points[:, 1],
            sampled,
        )

    def _update_selection(self):
        self.status_artist.set_text(
            self._status_message()
        )

        if self.start is None or self.end is None:
            self.fig.canvas.draw_idle()
            return

        geometry = self._selection_geometry()

        if geometry is None:
            self.centerline_artist.set_data(
                [self.start[0]],
                [self.start[1]],
            )
            self.centerline_artist.set_visible(True)
            self.fig.canvas.draw_idle()
            return

        length, tangent, normal, corners = geometry

        (
            self.distance,
            self.profile,
            self.x_line,
            self.y_line,
            self.sample_values,
        ) = self._extract_profile(
            length,
            tangent,
            normal,
        )

        # Update selection artists.
        self.centerline_artist.set_data(
            [self.start[0], self.end[0]],
            [self.start[1], self.end[1]],
        )
        self.centerline_artist.set_visible(True)

        self.strip_artist.set_xy(corners)
        self.strip_artist.set_visible(True)

        # Update profile plot.
        self.profile_artist.set_data(
            self.distance,
            self.profile,
        )

        self.profile_ax.relim()
        self.profile_ax.autoscale_view()
        self.profile_ax.set_title(
            f"Averaged lineout — width = {self.width:.4g}"
        )

        self.fig.canvas.draw_idle()
        self.profile_fig.canvas.draw_idle()

    def get_result(self, copy=True):
        """
        Return the latest lineout.

        Returns
        -------
        distance, profile, x_line, y_line
        """
        if self.distance is None:
            raise RuntimeError(
                "No lineout has been selected yet."
            )

        result = (
            self.distance,
            self.profile,
            self.x_line,
            self.y_line,
        )

        if copy:
            return tuple(array.copy() for array in result)

        return result

    def save(self, filename, delimiter=","):
        """
        Save distance, x, y, and averaged profile to a text file.
        """
        distance, profile, x_line, y_line = self.get_result()

        output = np.column_stack(
            [
                distance,
                x_line,
                y_line,
                profile,
            ]
        )

        header = delimiter.join(
            [
                "distance",
                "x",
                "y",
                "profile",
            ]
        )

        np.savetxt(
            filename,
            output,
            delimiter=delimiter,
            header=header,
            comments="",
        )

    def clear(self):
        self.start = None
        self.end = None
        self.distance = None
        self.profile = None
        self.x_line = None
        self.y_line = None
        self.sample_values = None

        self.centerline_artist.set_visible(False)
        self.strip_artist.set_visible(False)

        self.profile_artist.set_data([], [])
        self.profile_ax.relim()
        self.profile_ax.autoscale_view()

        self.fig.canvas.draw_idle()
        self.profile_fig.canvas.draw_idle()

    def disconnect(self):
        """
        Disconnect the mouse and keyboard callbacks.
        """
        for callback_id in self._callback_ids:
            self.fig.canvas.mpl_disconnect(callback_id)

        self._callback_ids.clear()