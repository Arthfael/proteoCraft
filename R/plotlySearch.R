#' .plotlySearch
#' 
#' @description
#' This function adds a search field to a plotly volcano plot.
#' It is NOT my work: it was designed under my supervision by chatGPT and istaGPT (our local GPT).
#' 
#' Note that this function has special behaviour for lightgrey points (rgba = 211, 211, 211, 1): matching points are set to rgba = 0, 0, 0, 1 while non matching points are set to rgba = 211, 211, 211, 0.2!
#' This allows seeing contrast between selected/non-selected points for the grey ("boring") area of volcano plots, where we use "lightgrey'.
#' 
#' @param p Plotly scatter plot (built or not) which should receive a search field.
#' @param search_data A data.frame derived from the data used to create the trace(s) of interest, with the same rows (number and ordering) and searchable columns.
#' @param fade_opacity Opacity of non-selected points. Default = 0.2.
#' @param placeholder Text displayed in the search field, default = "Search..."
#' 
#' @returns
#' A build plotly widget with the search field added.
#' 
#' @examples
#' plotLy <- do.call(plotly::plot_ly, args_ly_)
#' plotLy <- .plotlySearch(plotLy, searchDat)
#' 
#' @export


.plotlySearch <- function(p,
                          search_data,
                          fade_opacity = 0.2,
                          placeholder = "Search...") {
  #`%||%` <- \(x, y) { if (is.null(x) || (!length(x))) { y } else { x } }
  p <- plotly::plotly_build(p)
  # Identify scatter/scattergl marker traces
  is_marker_trace <- vapply(p$x$data, \(tr) {
    tr_type <- as.character(tr$type %||% "scatter")[1L]
    tr_mode <- as.character(tr$mode %||% "")[1L]
    has_xy <- (!is.null(tr$x)) && (!is.null(tr$y))
    is_scatter <- tr_type %in% c("scatter", "scattergl")
    has_markers <- grepl("markers", tr_mode, fixed = TRUE)
    return(has_xy && is_scatter && has_markers)
  }, TRUE)
  trace_indices <- which(is_marker_trace)
  if (!length(trace_indices)) {
    stop("No scatter/scattergl marker traces found.")
  }
  n_points <- vapply(trace_indices, \(i) { length(p$x$data[[i]]$x) }, 1L)
  total_points <- sum(n_points)
  search_data <- as.data.frame(search_data)
  if (nrow(search_data) != total_points) {
    warning("Search_data has ", nrow(search_data), " rows, but the selected marker traces contain ", total_points, " points in total -> no search bar added!")
    return(p)
  }
  # Ensure marker object exists
  for (i in trace_indices) {
    if (is.null(p$x$data[[i]]$marker)) {
      p$x$data[[i]]$marker <- list()
    }
  }
  search_json <- jsonlite::toJSON(search_data,
                                  dataframe = "rows",
                                  auto_unbox = TRUE,
                                  na = "null")
  trace_json <- jsonlite::toJSON(as.integer(trace_indices - 1L),
                                 auto_unbox = FALSE)
  n_points_json <- jsonlite::toJSON(as.integer(n_points),
                                    auto_unbox = FALSE)
  fade_json <- jsonlite::toJSON(fade_opacity,
                                auto_unbox = TRUE)
  placeholder_json <- jsonlite::toJSON(placeholder,
                                       auto_unbox = TRUE)
  js <- sprintf("function(el, x) {
  const gd = el;
  const searchData = %s;
  const traceIndices = %s;
  const nPoints = %s;
  const fadeOpacity = %s;
  const placeholder = %s;
  const greyColorReplacement = 'rgba(0, 0, 0, 1)';
  function isTargetGreyColor(value) {
    if (value === null || value === undefined || typeof value !== 'string') {
      return false;
    }
    const v = value.toLowerCase().replace(/\\s+/g, '');
    return (v === 'rgba(211,211,211,1)' || v === 'rgb(211,211,211,1)' || v === 'rgb(211,211,211)');
  }
  function replaceGreyWithWhite(color, n) {
    if (Array.isArray(color)) {
      return color.map(function(c) {
        return isTargetGreyColor(c)
          ? greyColorReplacement
          : c;
      });
    }
    if (color !== null && color !== undefined && typeof color.length === 'number' && typeof color !== 'string' && color.length === n) {
      return Array.prototype.slice.call(color).map(function(c) {
        return isTargetGreyColor(c)
          ? greyColorReplacement
          : c;
      });
    }
    if (isTargetGreyColor(color)) {
      return Array(n).fill(greyColorReplacement);
    }
    return color;
  }
  // Search UI
  const wrapper = document.createElement('div');
  const input = document.createElement('input');
  input.type = 'text';
  input.placeholder = placeholder;
  input.style.width = '250px';
  input.style.padding = '5px 8px';
  input.style.marginBottom = '8px';
  input.style.border = '1px solid #ccc';
  input.style.borderRadius = '4px';
  const count = document.createElement('span');
  count.style.marginLeft = '10px';
  wrapper.appendChild(input);
  wrapper.appendChild(count);
  el.parentNode.insertBefore(wrapper, el);
  // Store original marker opacity and color per trace
  const originalOpacityByTrace = traceIndices.map(function(traceIndex, j) {
    const tr = gd.data[traceIndex];
    const n = nPoints[j];
    if (!tr.marker) {
      tr.marker = {};
    }
    const op = tr.marker.opacity;
    if (Array.isArray(op)) {
      return op.slice();
    }
    if (op !== null && op !== undefined && typeof op.length === 'number' && typeof op !== 'string' && op.length === n) {
      return Array.prototype.slice.call(op);
    }
    if (op === null || op === undefined) {
      return null;
    }
    return op;
  });
  const originalColorByTrace = traceIndices.map(function(traceIndex, j) {
    const tr = gd.data[traceIndex];
    const n = nPoints[j];
    if (!tr.marker) {
      tr.marker = {};
    }
    const col = tr.marker.color;
    if (Array.isArray(col)) {
      return col.slice();
    }
    if (col !== null && col !== undefined && typeof col.length === 'number' && typeof col !== 'string' && col.length === n) {
      return Array.prototype.slice.call(col);
    }
    if (col === null || col === undefined) {
      return null;
    }
    return col;
  });
  input.addEventListener('keydown', function(event) {
    if (event.key !== 'Enter') {
      return;
    }
    const query = input.value.trim().toLowerCase();
    // Empty search: restore original colors and original opacity
    if (query === '') {
      traceIndices.forEach(function(traceIndex, j) {
        Plotly.restyle(gd, {
            'marker.opacity': [originalOpacityByTrace[j]],
            'marker.color': [originalColorByTrace[j]]
          }, [traceIndex]);
      });
      count.textContent = '';
      return;
    }
    let offset = 0;
    let totalMatches = 0;
    traceIndices.forEach(function(traceIndex, j) {
      const n = nPoints[j];
      const rowsForTrace = searchData.slice(offset, offset + n);
      offset += n;
      const matches = rowsForTrace.map(function(row) {
        return Object.values(row).some(function(value) {
          if (value === null || value === undefined) {
            return false;
          }
          return String(value).toLowerCase().includes(query);
        });
      });
      totalMatches += matches.filter(Boolean).length;
      // Absolute opacity:
      // matching points stay fully opaque,
      // non-matching points are faded.
      const newOpacity = matches.map(function(isMatch) {
        return isMatch ? 1 : fadeOpacity;
      });
      // Temporarily convert original light-grey points to white
      const newColor = replaceGreyWithWhite(originalColorByTrace[j], n);
      Plotly.restyle(gd, {
          'marker.color': [newColor],
          'marker.opacity': [newOpacity]
        }, [traceIndex]);
    });
    count.textContent = totalMatches + (totalMatches === 1 ? ' match' : ' matches');
  });
}",
                search_json,
                trace_json,
                n_points_json,
                fade_json,
                placeholder_json)
  return(htmlwidgets::onRender(p, js))
}
