using ColorSchemes, Colors, CairoMakie
using LaTeXStrings

function simple_theme!()
  set_theme!(theme_minimal(); figure_padding=16, fontsize=20)
  cmap = ColorScheme([RGBA(1, 1, 1, 0.0), RGBA(0.4, 0.8, 0.0, 0.5), RGBA(0, 0.5, 0, 0.5), RGBA(1, 1, 0, 0.5)])

  theme = Theme(;
    Figure=(backgroundcolor = :transparent),
    Heatmap=(colormap=amp, lowclip=:transparent),
    Scatter=(markersize=20, colormap=cmap),
  )
  update_theme!(theme)
  return theme
end

@deprecate customTheme!(fs=28) default_theme!(fs=28)
customTheme!(fs=28) = default_theme!(fs)

function scalebar!(ax, position; width=100, height=10, fontsize=28, color=:white, offset=20)
  x, y = position
  lines!(ax, [x - width / 2, x + width / 2], [y - height, y - height]; linewidth=10, color=color)
  return text!(
    ax, L"%$width"; position=position .+ (0, offset), align=(:center, :center), fontsize=fontsize, color=color
  )
end

function default_theme!(fs=28)
  set_theme!(theme_minimal(); fontsize=fs)

  theme = Theme(;
    Axis=(
      xgridvisible=false,
      ygridvisible=false,
      xtrimspine=false,
      ytrimspine=(false, false),
      xlabelsize=40,
      ylabelsize=50,
      backgroundcolor=:transparent,
    ),
  )
  update_theme!(theme)
  return theme
end

rect(x, y, width, aspect) = (x, x + width, y, y + width * aspect)

function pointPlot!(ax, x, y, x_err, y_err, config)
  scatter!(ax, [x], [y]; config...)
  errorbars!(ax, [x], [y], [y_err]; whiskerwidth=12, color=config.color, direction=:y)
  errorbars!(ax, [x], [y], [x_err]; whiskerwidth=12, color=config.color, direction=:x)
  return nothing
end

function alphaColor(cscheme::ColorScheme, alpha=0.5; ncolors=12)
  return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1; length=ncolors)])
end

function addtext!(ax, txt, pos, rot, fs=25)
  return text!(ax, txt; position=pos, align=(:left, :center), rotation=rot, fontsize=fs)
end

function addline!(ax, vlow, vhigh; shift=1, pow=0, color=:black)
  vrange = LinRange(vlow, vhigh, 100)
  return lines!(ax, vrange, (10.0 .^ shift) .* vrange .^ pow; color=color, linewidth=3)
end