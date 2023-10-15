
from django.contrib import admin
from django.urls import path
from prvi.views import Index
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from . import views



urlpatterns=[
    path('admin/', admin.site.urls),
    path('', Index.as_view(), name='index'),
    path('preuzimanje/', views.Preuzimanje, name='preuzimanje'),
    
]


urlpatterns += staticfiles_urlpatterns()
