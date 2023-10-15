

from django import forms

class UploadForm(forms.Form):
    file = forms.FileField()

class XYZForm(forms.Form):
    x1 = forms.IntegerField(label='X1')
    y1 = forms.IntegerField(label='Y1')
    x2 = forms.IntegerField(label='X2')
    y2 = forms.IntegerField(label='Y2')
    x3 = forms.IntegerField(label='X3')
    y3 = forms.IntegerField(label='Y3')
    x4 = forms.IntegerField(label='X4')
    y4 = forms.IntegerField(label='Y4')

class BCForm(forms.Form):
    BCx1 = forms.IntegerField(label='BCX1')
    BCy1 = forms.IntegerField(label='BCY1')
    BCx2 = forms.IntegerField(label='BCX2')
    BCy2 = forms.IntegerField(label='BCY2')
    BCx3 = forms.IntegerField(label='BCX3')
    BCy3 = forms.IntegerField(label='BCY3')
    BCx4 = forms.IntegerField(label='BCX4')
    BCy4 = forms.IntegerField(label='BCY4')

class FORCEForm(forms.Form):
    Fx1 = forms.IntegerField(label='FX1')
    Fy1 = forms.IntegerField(label='FY1')
    Fx2 = forms.IntegerField(label='FX2')
    Fy2 = forms.IntegerField(label='FY2')
    Fx3 = forms.IntegerField(label='FX3')
    Fy3 = forms.IntegerField(label='FY3')
    Fx4 = forms.IntegerField(label='FX4')
    Fy4 = forms.IntegerField(label='FY4')

class ELDATAForm(forms.Form):
    YE = forms.FloatField(label='YE')
    nu = forms.FloatField(label='nu')
    Debljina = forms.FloatField(label='Debljina')

