---
layout: page
title: "Download"
group: navigation
---

{% include JB/setup %}

If using a Mac, the easiest way to download and install __kallisto__ is via brew with the commands:

```
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install kallisto
```

Alternatively, with [conda](https://conda.io/docs/index.html) installed and the [bioconda](https://bioconda.github.io/index.html#set-up-channels) channel activated, __kallisto__ can be installed with:

```
conda install kallisto
```

The __kallisto__ bioconda installation will work with 64 bit linux or Mac OS. More details are available at the [__kallisto__ bioconda page](https://bioconda.github.io/recipes/kallisto/README.html).

__kllisto__ can also be installed on [FreeBSD](https://www.freebsd.org/ports/index.html) via the FreeBSD ports system using

```
pkg install kallisto
```

__kallisto__ binaries for Mac OS X, NetBSD, RHEL/CentOS and SmartOS can be installed on most POSIX platforms using [pkgsrc](https://www.pkgsrc.org):

```
pkgin install kallisto
```

For other operating systems download executables or source from the links below and follow the provided instructions. Note that prebuilt kallisto indices from the human transcriptome and many model organism transcriptomes are available from the [kallisto transcriptome indices](https://github.com/pachterlab/kallisto-transcriptome-indices) page.

#### Releases

The __kallisto__ GitHub repository is [here](http://github.com/pachterlab/kallisto).

<table class="table">
  <thead>
    <tr>
      <th style="text-align: left">Version</th>
      <th>Date</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>

{% for post in site.categories.releases %}
    <tr>
    	<td>Release notes: <a href="{{ site.url }}/kallisto/{{ post.url }}">{{ post.version }}</a></td>
    	<td><span class="entry-date"><time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%B %d, %Y" }}</time></span></td>

        <td><a href="https://github.com/pachterlab/kallisto/releases/download/{{ post.version }}/kallisto_mac-{{ post.version }}.tar.gz">Mac</a></td>
        <td>{% if post.version != 'v0.42.2' %} <a href="https://github.com/pachterlab/kallisto/releases/download/{{ post.version }}/kallisto_linux-{{ post.version }}.tar.gz">Linux</a> {% endif %} </td>
        <td>{% if post.version >= 'v0.43.0' %} <a href="https://github.com/pachterlab/kallisto/releases/download/{{ post.version }}/kallisto_windows-{{ post.version }}.zip">Windows</a> {% endif %} </td>
        <td>{% if post.version >= 'v0.44.0' %} <a href="https://github.com/pachterlab/kallisto/releases/download/{{ post.version }}/kallisto_rock64-{{ post.version }}.tar.gz">Rock64</a> {% endif %} </td>
        <td><a href="https://github.com/pachterlab/kallisto/archive/{{ post.version }}.tar.gz">Source</a></td>
    </tr>
{% endfor %}
</table>

#### Licence

kallisto is distributed under the BSD 2-clause "Simplified" License:

~~~
BSD 2-Clause License

Copyright (c) 2017, Nicolas Bray, Harold Pimentel, Páll Melsted and Lior Pachter
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
~~~
