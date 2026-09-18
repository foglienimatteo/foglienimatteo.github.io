Posters and other PDFs shown as a preview card on /Publications/.

To add one:

1. drop the PDF here, e.g. `2025-06-pasc25-poster.pdf`
2. `bin/pdf-preview assets/posters/2025-06-pasc25-poster.pdf`
   (writes `2025-06-pasc25-poster.jpg` next to it)
3. reference both from `_pages/Work.md`:

       {% include pdf-embed.html
            pdf="/assets/posters/2025-06-pasc25-poster.pdf"
            image="/assets/posters/2025-06-pasc25-poster.jpg"
            alt="..."
            caption="..." %}

The include renders nothing until the PDF actually exists, so an unfinished
entry never shows a broken image on the live site.

Poster PDFs exported from Illustrator or Inkscape are often 50-200MB. Git
keeps every binary forever, so shrink it first, e.g.

    gs -sDEVICE=pdfwrite -dPDFSETTINGS=/ebook -dNOPAUSE -dQUIET -dBATCH \
       -sOutputFile=small.pdf big.pdf
