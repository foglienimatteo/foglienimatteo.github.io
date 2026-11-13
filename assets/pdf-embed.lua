--[[
  Clickable preview of a PDF: a raster of page 1 that opens the real PDF full
  page in a new tab.

      {{< pdf-embed pdf="/assets/posters/ISC26_Poster.pdf"
                    image="/assets/posters/ISC26_Poster.jpg"
                    alt="ISC26 poster: Rethinking Compute Platforms"
                    caption="ISC High Performance 2026"
                    width="36rem" >}}

  Generate `image` from `pdf` with:  bin/pdf-preview <path-to-pdf>

  Why a static image and not an <object>/<iframe> of the PDF itself: iOS
  Safari and Android Chrome will not render a PDF inline, and an embedded PDF
  viewer swallows pointer events, so wrapping it in a link does not give
  click-to-open. A JPEG is a fraction of the size and renders everywhere.

  Paths start with / and are relative to the site root.
]]

local function esc(s)
  return (s:gsub("&", "&amp;"):gsub("<", "&lt;"):gsub(">", "&gt;"):gsub('"', "&quot;"))
end

local function arg(kwargs, name)
  local v = kwargs[name]
  if v == nil then return "" end
  return pandoc.utils.stringify(v)
end

return {
  ["pdf-embed"] = function(args, kwargs)
    local pdf, image = arg(kwargs, "pdf"), arg(kwargs, "image")
    local alt, caption = arg(kwargs, "alt"), arg(kwargs, "caption")
    local width, height = arg(kwargs, "width"), arg(kwargs, "height")

    local box = ""
    if width ~= "" then box = box .. "--pdf-w: " .. width .. ";" end
    if height ~= "" then
      box = box .. "--pdf-h: " .. height .. ";"
    elseif width ~= "" then
      box = box .. "--pdf-h: none;"
    end
    local style = box ~= "" and (' style="' .. esc(box) .. '"') or ""

    local label = alt ~= "" and alt or "Open the PDF"
    local cap = caption ~= "" and (esc(caption) .. " &middot; ") or ""

    local html = string.format([[
<figure class="pdf-embed"%s>
  <a class="pdf-embed__preview" href="%s" target="_blank" rel="noopener" aria-label="%s (opens full page)">
    <img src="%s" alt="%s" loading="lazy" decoding="async">
  </a>
  <figcaption>%s<a href="%s" target="_blank" rel="noopener">open full page (PDF)</a></figcaption>
</figure>]], style, esc(pdf), esc(label), esc(image), esc(alt), cap, esc(pdf))

    return pandoc.RawBlock("html", html)
  end
}
