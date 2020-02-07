from django.shortcuts import render

posts = [
    {
        'author': 'CoreyMS',
        'title': 'Blog Post 1',
        'content': 'First post content',
        'date_posted': 'August 27, 2018'
    },
    {
        'author': 'Jane Doe',
        'title': 'Blog Post 2',
        'content': 'Second post content',
        'date_posted': 'August 28, 2018'
    }
]


def startpage(request):
    context = {'posts': posts}
    return render(request, 'startpage/home.html', context)


def about(request):
    return render(request, 'startpage/about.html', context={'title': 'about autometa'})
