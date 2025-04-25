---
title: "All options of posts"
#date:
#excerpt_separator: "<!--more-->"
last_modified_at: 2025-04-25T16:20:02-05:00
categories:
  - Blog
tags:
  - Post Formats
  - readability
  - standard
# link: https://github.com
---



# All options of Posts 

- [All options of Posts](#all-options-of-posts)
    - [Introduction](#introduction)
    - [last\_modified\_at](#last_modified_at)
    - [excerpt\_separator](#excerpt_separator)
    - [Link](#link)
    - [Quotes](#quotes)
    - [Notice](#notice)
      - [New Site Features](#new-site-features)

<br>
<br>
<br>



### Introduction

You'll find this post in your `_posts` directory. Go ahead and edit it and re-build the site to see your changes. You can rebuild the site in many different ways, but the most common way is to run `jekyll serve`, which launches a web server and auto-regenerates your site when a file is updated.

To add new posts, simply add a file in the `_posts` directory that follows the convention `YYYY-MM-DD-name-of-post.ext` and includes the necessary front matter. Take a look at the source for this post to get an idea about how it works.

Jekyll also offers powerful support for code snippets:

```ruby
def print_hi(name)
  puts "Hi, #{name}"
end
print_hi('Tom')
#=> prints 'Hi, Tom' to STDOUT.
```

Check out the [Jekyll docs][jekyll-docs] for more info on how to get the most out of Jekyll. File all bugs/feature requests at [Jekyll’s GitHub repo][jekyll-gh]. If you have questions, you can ask them on [Jekyll Talk][jekyll-talk].

[jekyll-docs]: https://jekyllrb.com/docs/home
[jekyll-gh]:   https://github.com/jekyll/jekyll
[jekyll-talk]: https://talk.jekyllrb.com/

<br>
<br>
<br>







### last_modified_at

This post has been updated and should show a modified date if used in a layout.

Mrs. Darling first heard of Peter when she was tidying up her children's minds. It is the nightly custom of every good mother after her children are asleep to rummage in their minds and put things straight for next morning, repacking into their proper places the many articles that have wandered during the day.

<!--more-->

<br>
<br>
<br>







### excerpt_separator

This post has a manual excerpt `<!--more-->` set after the second paragraph. The following YAML Front Matter has also be applied:

```yaml
excerpt_separator: "<!--more-->"
```


<br>
<br>
<br>




### Link


This theme supports **link posts**, made famous by John Gruber. To use, just add `link: http://url-you-want-linked` to the post's YAML front matter and you're done.

> And this is how a quote looks.

Some [link](#) can also be shown.



<br>
<br>
<br>




### Quotes


> Only one thing is impossible for God: To find any sense in any copyright law on the planet.
  
> <cite><a href="http://www.brainyquote.com/quotes/quotes/m/marktwain163473.html">Mark Twain</a></cite>


<br>
<br>
<br>






### Notice


A notice displays information that explains nearby content. Often used to call attention to a particular detail.

When using Kramdown `{: .notice}` can be added after a sentence to assign the `.notice` to the `<p></p>` element. 

**Changes in Service:** We just updated our [privacy policy](#) here to better service our customers. We recommend reviewing the changes.
{: .notice}

**Primary Notice:** Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer nec odio. [Praesent libero](#).  
{: .notice--primary}

**Info Notice:** Lorem ipsum dolor sit amet, [consectetur adipiscing elit](#). 
{: .notice--info}

**Warning Notice:** Lorem ipsum dolor sit amet, consectetur adipiscing elit. [Integer nec odio](#). Praesent libero.  
{: .notice--warning}

**Danger Notice:** Lorem ipsum dolor sit amet, [consectetur adipiscing](#) elit. 
{: .notice--danger}

**Success Notice:** Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nulla quis sem at [nibh elementum](#) imperdiet.
{: .notice--success}

Want to wrap several paragraphs or other elements in a notice? Using Liquid to capture the content and then filter it with `markdownify` is a good way to go.

```html
{% raw %}{% capture notice-2 %}
#### New Site Features

* You can now have cover images on blog pages
* Drafts will now auto-save while writing
{% endcapture %}{% endraw %}

<div class="notice">{% raw %}{{ notice-2 | markdownify }}{% endraw %}</div>
```

{% capture notice-2 %}
#### New Site Features

* You can now have cover images on blog pages
* Drafts will now auto-save while writing
{% endcapture %}

<div class="notice">
  {{ notice-2 | markdownify }}
</div>

Or you could skip the capture and stick with straight HTML.

```html
<div class="notice">
  <h4>Message</h4>
  <p>A basic message.</p>
</div>
```

<div class="notice">
  <h4>Message</h4>
  <p>A basic message.</p>
</div>