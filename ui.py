# -*- coding: utf-8 -*-
import wx
import wx.grid

import matplotlib
matplotlib.interactive(False)
matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

import numpy as np

import json


G = 6.67384


def to_str(value):
    if value is None:
        return ''
    else:
        return str(value)


def to_int(value):
    if value:
        return int(float(value))
    else:
        return None


def to_float(value):
    if value:
        return float(value)
    else:
        return None


class Prism(object):
    def __init__(self, rho=None, lx=None, ly=None, zt=None, zb=None, xc=None, yc=None, phi=None):
        self.rho = rho
        self.lx = lx
        self.ly = ly
        self.zt = zt
        self.zb = zb
        self.xc = xc
        self.yc = yc
        self.phi = phi


class Grid(object):
    def __init__(self, nx=None, ny=None, dx=None, dy=None, xo=None, yo=None):
        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy
        self.xo = xo
        self.yo = yo


def readfile(filepath):
    gridinfotransforms = [int, int, float, float, float, float]
    with open(filepath, 'r') as f:
        
        gridinfo = [transform(value) for transform, value in zip(gridinfotransforms, f.readline().split())]
        bodyinfo = [[float(b) for b in a.split()] for a in f.readlines()]
    
    return gridinfo, bodyinfo


def writefile(filepath, grid, prisms):
    gridformat = '{nx:g} {ny:g} {dx:g} {dy:g} {xo:g} {yo:g}'
    prismformat = '{rho:g} {lx:g} {ly:g} {zt:g} {zb:g} {xc:g} {yc:g} {phi:g}'
    
    with open(filepath, 'w') as f:
        print >> f, gridformat.format(**grid.__dict__)
        
        for prism in prisms:
            print >> f, prismformat.format(**prism.__dict__)


def gkernel(x, y, z):
    r = np.sqrt(x*x + y*y + z*z)

    if z != 0.0:
        return y*np.log(r + x) + x*np.log(r + y) + 2.0*z*np.arctan((x + y + r)/z)

    k = np.zeros_like(x)

    bx = (x != 0.0)
    by = (y != 0.0)

    b1 = bx*(~by)
    b2 = (~bx)*by
    b3 = bx*by

    if b1.any():
        k[b1] = x[b1]*np.log(r[b1])

    if b2.any():
        k[b2] = y[b2]*np.log(r[b2])

    if b3.any():
        k[b3] = y[b3]*np.log(r[b3] + x[b3]) + \
                x[b3]*np.log(r[b3] + y[b3])

    return k


def gfield(X, Y, prism):
    cos_phi = np.cos(prism.phi*np.pi/180.0)
    sin_phi = np.sin(prism.phi*np.pi/180.0)
    dx = (prism.xc - X)*cos_phi - (prism.yc - Y)*sin_phi
    dy = (prism.xc - X)*sin_phi + (prism.yc - Y)*cos_phi
    x1 = dx - prism.lx
    x2 = dx + prism.lx
    y1 = dy - prism.ly
    y2 = dy + prism.ly
    z1 = prism.zt
    z2 = prism.zb
    k = gkernel(x1, y1, z1) - gkernel(x2, y1, z1) - gkernel(x1, y2, z1) + gkernel(x2, y2, z1)
    if z1 != z2:
        k -= gkernel(x1, y1, z2) - gkernel(x2, y1, z2) - gkernel(x1, y2, z2) + gkernel(x2, y2, z2)
    return G*prism.rho*k


def batchgfield(X, Y, prisms):
    nx, ny = X.shape
    nprisms = len(prisms)
    field = np.zeros((nx, ny, nprisms), dtype=float)
    for i in range(nprisms):
        field[:,:,i] = gfield(X, Y, prisms[i])
    return np.sum(field, axis=2)


# class SaveFigDialog(wx.Dialog):
    # default_props = {'h': 6.0, 'l': 6.0, 'dpi': 100.0}
    
    # def __init__(self, *args, **kwargs):
        # super(SaveFigDialog, self).__init__(*args, **kwargs)
        
        # fig_fgs = wx.FlexGridSizer(2, 2, 5, 5)
        
        # size_label = wx.StaticText(self, label=u"Tamanho (pol):")
        # dpi_label = wx.StaticText(self, label=u"Resolução (dpi):")
        
        # self.h_ctrl = wx.TextCtrl(self)
        # self.l_ctrl = wx.TextCtrl(self)
        # self.dpi_ctrl = wx.TextCtrl(self)
        
        # size_sizer = wx.BoxSizer(wx.HORIZONTAL)
        # size_sizer.Add(self.h_ctrl, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, 5)
        # size_sizer.Add(self.l_ctrl, 1, wx.EXPAND | wx.RIGHT, 5)
        
        # fig_fgs.AddMany([(size_label), (size_sizer, 1, wx.EXPAND),
                         # (dpi_label), (self.dpi_ctrl)])
        
        # fig_fgs.AddGrowableCol(1, 1)
        
        # button_sizer = self.CreateButtonSizer(wx.OK | wx.CANCEL)

        # sizer = wx.BoxSizer(wx.VERTICAL)
        # sizer.Add(fig_fgs, 1, wx.EXPAND | wx.ALL, 5)
        # sizer.Add(button_sizer, 0, wx.EXPAND | wx.ALL, 5)

        # self.SetSizer(sizer)
        # self.Fit()
        # self.SetTitle(u"Opções da figura")
        
        # self.set_props(SaveFigDialog.default_props.copy())
    
    # def set_props(self, props):
        # h = props.get('h', None)
        # l = props.get('l', None)
        # dpi = props.get('dpi', None)
        
        # self.h_ctrl.SetValue(to_str(h))
        # self.l_ctrl.SetValue(to_str(l))
        # self.dpi_ctrl.SetValue(to_str(dpi))
    
    # def get_props(self):
        # props = {}
        
        # h = self.h_ctrl.GetValue().strip()
        # l = self.l_ctrl.GetValue().strip()
        # dpi = self.dpi_ctrl.GetValue().strip()
        
        # props['h'] = to_float(h)
        # props['l'] = to_float(l)
        # props['dpi'] = to_float(dpi)
        
        # return props


class ConfigDialog(wx.Dialog):
    # cmaps_pt =  [u'Cinza', u'Arco-íris', u'Jet', u'HSV', u'Frio', u'Quente', u'Cobre', u'Osso', u'Rosa', u'Oceano', u'Primavera', u'Verão', u'Outono', u'Inverno']
    cmaps_pt = ['Gray', 'Rainbow', 'Jet', 'HSV', 'Cool', 'Hot', 'Copper', 'Bone', 'Pink', 'Ocean', 'Spring', 'Summer', 'Autumn', 'Winter']
    cmaps = ['gray', 'rainbow', 'jet', 'hsv', 'cool', 'hot', 'copper', 'bone', 'pink', 'ocean', 'spring', 'summer', 'autumn', 'winter']
    colors_pt = [u'Preto', u'Branco', u'Vermelho', u'Verde', u'Azul', u'Ciano', u'Magenta', u'Amarelo']
    colors = ['k', 'w', 'r', 'g', 'b', 'c', 'm', 'y']
    
    def __init__(self, *args, **kwargs):
        super(ConfigDialog, self).__init__(*args, **kwargs)

        ### GRID ###
        grid_sb = wx.StaticBox(self, label=u"Grid:")
        grid_sbs = wx.StaticBoxSizer(grid_sb, wx.HORIZONTAL)

        grid_fgs = wx.FlexGridSizer(2, 6, 5, 5)
        nx_label = wx.StaticText(self, label="nx:")
        dx_label = wx.StaticText(self, label="dx:")
        xo_label = wx.StaticText(self, label="xo:")
        ny_label = wx.StaticText(self, label="ny:")
        dy_label = wx.StaticText(self, label="dy:")
        yo_label = wx.StaticText(self, label="yo:")
        self.nx_ctrl = wx.TextCtrl(self)
        self.dx_ctrl = wx.TextCtrl(self)
        self.xo_ctrl = wx.TextCtrl(self)
        self.ny_ctrl = wx.TextCtrl(self)
        self.dy_ctrl = wx.TextCtrl(self)
        self.yo_ctrl = wx.TextCtrl(self)

        grid_fgs.AddMany([(nx_label), (self.nx_ctrl, 1, wx.EXPAND), (dx_label), (self.dx_ctrl, 1, wx.EXPAND), (xo_label), (self.xo_ctrl, 1, wx.EXPAND),
                          (ny_label), (self.ny_ctrl, 1, wx.EXPAND), (dy_label), (self.dy_ctrl, 1, wx.EXPAND), (yo_label), (self.yo_ctrl, 1, wx.EXPAND)])

        grid_fgs.AddGrowableCol(1, 1)
        grid_fgs.AddGrowableCol(3, 1)
        grid_fgs.AddGrowableCol(5, 1)

        grid_sbs.Add(grid_fgs, 1, wx.EXPAND)
        ### GRID ###
        ### FIG ###
        fig_sb = wx.StaticBox(self, label=u"Figura:")
        fig_sbs = wx.StaticBoxSizer(fig_sb, wx.HORIZONTAL)

        fig_fgs = wx.FlexGridSizer(8, 4, 5, 5)

        title_label = wx.StaticText(self, label=u"Título:")
        xlabel_label = wx.StaticText(self, label=u"Legenda x:")
        ylabel_label = wx.StaticText(self, label=u"Legenda y:")
        pos_label = wx.StaticText(self, label=u"Posição:")
        contour_label = wx.StaticText(self, label=u"Contorno:")
        color_label = wx.StaticText(self, label=u"Cor:")
        levels_label = wx.StaticText(self, label=u"Curvas:")
        clabel_label = wx.StaticText(self, label=u"Valores nas curvas:")
        fill_label = wx.StaticText(self, label=u"Imagem:")
        cmap_label = wx.StaticText(self, label=u"Escala de cor:")
        cbar_label = wx.StaticText(self, label=u"Barra de cor:")

        self.title_ctrl = wx.TextCtrl(self)
        self.xlabel_ctrl = wx.TextCtrl(self)
        self.ylabel_ctrl = wx.TextCtrl(self)
        pos_szr = wx.BoxSizer(wx.HORIZONTAL)
        self.pos_ctrls = []
        for i in range(4):
            pos_ctrl = wx.TextCtrl(self, size=(25, -1))
            self.pos_ctrls.append(pos_ctrl)
            pos_szr.Add(pos_ctrl, 1, wx.EXPAND)
        self.contour_ctrl = wx.CheckBox(self)
        self.contour_ctrl.Bind(wx.EVT_CHECKBOX, self.on_contour)
        self.color_ctrl = wx.Choice(self, size=(25, -1))
        self.color_ctrl.AppendItems(self.colors_pt)
        self.levels_ctrl = wx.TextCtrl(self)
        self.clabel_ctrl = wx.CheckBox(self)
        self.fill_ctrl = wx.CheckBox(self)
        self.fill_ctrl.Bind(wx.EVT_CHECKBOX, self.on_fill)
        self.cmap_ctrl = wx.Choice(self, size=(25, -1))
        self.cmap_ctrl.AppendItems(self.cmaps_pt)
        self.cbar_ctrl = wx.CheckBox(self)

        label_style = wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT

        fig_fgs.AddMany([(title_label, 0, label_style), (self.title_ctrl, 1, wx.EXPAND), (0, 0), (0, 0),
                         (xlabel_label, 0, label_style), (self.xlabel_ctrl, 1, wx.EXPAND), (ylabel_label, 0, label_style), (self.ylabel_ctrl, 1, wx.EXPAND),
                         (pos_label, 0, label_style), (pos_szr, 1, wx.EXPAND), (0, 0), (0, 0),
                         (contour_label, 0, label_style), (self.contour_ctrl, 1, wx.EXPAND), (color_label, 0, label_style), (self.color_ctrl, 1, wx.EXPAND),
                         (0, 0), (0, 0), (levels_label, 0, label_style), (self.levels_ctrl, 1, wx.EXPAND),
                         (0, 0), (0, 0), (clabel_label, 0, label_style), (self.clabel_ctrl, 1, wx.EXPAND),
                         (fill_label, 0, label_style), (self.fill_ctrl, 1, wx.EXPAND), (cmap_label, 0, label_style), (self.cmap_ctrl, 1, wx.EXPAND),
                         (0, 0), (0, 0), (cbar_label, 0, label_style), (self.cbar_ctrl, 1, wx.EXPAND)])

        fig_fgs.AddGrowableCol(1, 1)
        fig_fgs.AddGrowableCol(3, 1)

        fig_sbs.Add(fig_fgs, 1, wx.EXPAND)
        ### FIG ###

        button_sizer = self.CreateButtonSizer(wx.OK | wx.CANCEL)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(grid_sbs, 0, wx.EXPAND)
        sizer.Add(fig_sbs, 0, wx.EXPAND)
        sizer.Add(button_sizer, 0, flag=wx.EXPAND)

        self.SetSizer(sizer)
        self.Fit()
        self.SetTitle(u"Configurações")

        self.on_contour(None)
        self.on_fill(None)

    def on_contour(self, event):
        self.color_ctrl.Enable(self.contour_ctrl.IsChecked())
        self.levels_ctrl.Enable(self.contour_ctrl.IsChecked())
        self.clabel_ctrl.Enable(self.contour_ctrl.IsChecked())

    def on_fill(self, event):
        self.cmap_ctrl.Enable(self.fill_ctrl.IsChecked())
        self.cbar_ctrl.Enable(self.fill_ctrl.IsChecked())

    def set_grid(self, grid):
        self.nx_ctrl.SetValue(to_str(grid.nx))
        self.dx_ctrl.SetValue(to_str(grid.dx))
        self.xo_ctrl.SetValue(to_str(grid.xo))
        self.ny_ctrl.SetValue(to_str(grid.ny))
        self.dy_ctrl.SetValue(to_str(grid.dy))
        self.yo_ctrl.SetValue(to_str(grid.yo))

    def get_grid(self):
        nx = self.nx_ctrl.GetValue()
        dx = self.dx_ctrl.GetValue()
        xo = self.xo_ctrl.GetValue()
        ny = self.ny_ctrl.GetValue()
        dy = self.dy_ctrl.GetValue()
        yo = self.yo_ctrl.GetValue()
        return Grid(to_int(nx), to_int(ny), to_float(dx), to_float(dy), to_float(xo), to_float(yo))

    def set_fig_props(self, props):
        self.title_ctrl.SetValue(props['title'])
        self.xlabel_ctrl.SetValue(props['xlabel'])
        self.ylabel_ctrl.SetValue(props['ylabel'])
        for i, pos in enumerate(props['pos']):
            self.pos_ctrls[i].SetValue(to_str(pos))
        self.contour_ctrl.SetValue(props['contour'])
        self.color_ctrl.SetSelection(self.colors.index(props['color']))
        if type(props['levels']) == list:
            self.levels_ctrl.SetValue(', '.join([str(a) for a in props['levels']]))
        else:
            self.levels_ctrl.SetValue(to_str(props['levels']))
        self.clabel_ctrl.SetValue(props['clabel'])
        self.fill_ctrl.SetValue(props['fill'])
        self.cmap_ctrl.SetSelection(self.cmaps.index(props['cmap']))
        self.cbar_ctrl.SetValue(props['cbar'])
        
        self.on_contour(None)
        self.on_fill(None)

    def get_fig_props(self):
        props = {}
        props['title'] = self.title_ctrl.GetValue()
        props['xlabel'] = self.xlabel_ctrl.GetValue()
        props['ylabel'] = self.ylabel_ctrl.GetValue()
        props['pos'] = [to_int(pos_ctrl.GetValue()) for pos_ctrl in self.pos_ctrls]
        props['contour'] = self.contour_ctrl.GetValue()
        props['color'] = self.colors[self.color_ctrl.GetSelection()]
        levels = self.levels_ctrl.GetValue()
        if ',' in levels:
            props['levels'] = [float(a) for a in levels.split(',')]
        else:
            props['levels'] = to_int(levels)
        props['clabel'] = self.clabel_ctrl.GetValue()
        props['fill'] = self.fill_ctrl.GetValue()
        props['cmap'] = self.cmaps[self.cmap_ctrl.GetSelection()]
        props['cbar'] = self.cbar_ctrl.GetValue()

        return props


class PrismTable(wx.grid.PyGridTableBase):
    properties = ["rho", "lx", "ly", "zt", "zb", "xc", "yc", "phi"]

    def __init__(self, prisms=[]):
        super(PrismTable, self).__init__()
        self.prisms = prisms

    def AppendRows(self, numRows=1):
        for i in range(numRows):
            self.prisms.append(Prism())

        self.GetView().BeginBatch()
        msg = wx.grid.GridTableMessage(self, wx.grid.GRIDTABLE_NOTIFY_ROWS_APPENDED, numRows)
        self.GetView().ProcessTableMessage(msg)
        self.GetView().EndBatch()

        return True

    def DeleteRows(self, pos=0, numRows=1):
        for i in range(pos, pos+numRows):
            prism = self.prisms.pop(i)
            del prism

        self.GetView().BeginBatch()
        msg = wx.grid.GridTableMessage(self, wx.grid.GRIDTABLE_NOTIFY_ROWS_DELETED, pos, numRows)
        self.GetView().ProcessTableMessage(msg)
        self.GetView().EndBatch()

        return True

    def GetColLabelValue(self, col):
        return PrismTable.properties[col]

    def GetRowLabelValue(self, row):
        return str(row + 1)

    def GetNumberCols(self):
        return 8

    def GetNumberRows(self):
        return len(self.prisms)

    def GetValue(self, row, col):
        value = getattr(self.prisms[row], PrismTable.properties[col], '')
        if value is None:
            value = ''
        return str(value)

    def SetValue(self, row, col, value):
        if value.strip():
            setattr(self.prisms[row], PrismTable.properties[col], float(value))
        else:
            setattr(self.prisms[row], PrismTable.properties[col], None)


class MainWindow(wx.Frame):
    iconsize = 32
    
    _figpropspath = 'figprops.json'

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.prisms = []
        self.xygrid = Grid()
        
        toolbar = self.CreateToolBar()
        toolbar.AddLabelTool(1, u'Abrir', wx.Bitmap('icons/open{}.png'.format(MainWindow.iconsize)), shortHelp=u'Abrir')
        toolbar.AddLabelTool(2, u'Salvar', wx.Bitmap('icons/save{}.png'.format(MainWindow.iconsize)), shortHelp=u'Salvar')
        toolbar.AddLabelTool(3, u'Salvar figura', wx.Bitmap('icons/savefig{}.png'.format(MainWindow.iconsize)), shortHelp=u'Salvar figura')
        toolbar.AddLabelTool(4, u'Exportar', wx.Bitmap('icons/export{}.png'.format(MainWindow.iconsize)), shortHelp=u'Exportar dados xyz')
        toolbar.AddLabelTool(5, u'Adicionar prisma', wx.Bitmap('icons/add{}.png'.format(MainWindow.iconsize)), shortHelp=u'Adicionar nova linha a tabela de prismas')
        toolbar.AddLabelTool(6, u'Remover prisma', wx.Bitmap('icons/remove{}.png'.format(MainWindow.iconsize)), shortHelp=u'Remover prismas selecionados')
        toolbar.AddLabelTool(7, u'Redesenhar', wx.Bitmap('icons/redraw{}.png'.format(MainWindow.iconsize)), shortHelp=u'Redesenhar a figura')
        toolbar.AddLabelTool(8, u'Configurações', wx.Bitmap('icons/config{}.png'.format(MainWindow.iconsize)), shortHelp=u'Configurações da figura')
        toolbar.Realize()

        self.Bind(wx.EVT_TOOL, self.on_open, id=1)
        self.Bind(wx.EVT_TOOL, self.on_save, id=2)
        self.Bind(wx.EVT_TOOL, self.on_savefig, id=3)
        self.Bind(wx.EVT_TOOL, self.on_export, id=4)
        self.Bind(wx.EVT_TOOL, self.on_add, id=5)
        self.Bind(wx.EVT_TOOL, self.on_remove, id=6)
        self.Bind(wx.EVT_TOOL, self.on_redraw, id=7)
        self.Bind(wx.EVT_TOOL, self.on_config, id=8)

        self.splitter = wx.SplitterWindow(self)

        self.table = PrismTable(self.prisms)

        self.grid = wx.grid.Grid(self.splitter)
        self.grid.SetDefaultColSize(60)
        self.grid.SetDefaultRowSize(20)
        self.grid.SetColLabelSize(20)
        self.grid.SetRowLabelSize(40)
        self.grid.SetTable(self.table)
        
        minsize = self.grid.GetEffectiveMinSize()[0] + 5

        self.figure = Figure()
        self.canvas = FigureCanvas(self.splitter, -1, self.figure)

        self.splitter.SplitVertically(self.canvas, self.grid)
        self.splitter.SetSashGravity(1.0)
        self.splitter.SetMinimumPaneSize(minsize)

        self.Maximize()
        
        sashpos = self.GetSize()[0] - minsize
        
        self.splitter.SetSashPosition(sashpos)
        
        self.SetTitle("GPrism")
        self.SetIcon(wx.Icon('icons/icon64.ico', wx.BITMAP_TYPE_ICO))

    def on_open(self, event):
        wildcard = "GPrism (*.gprism)|*.gprism"
        dlg = wx.FileDialog(self, "Abrir arquivo", "", "", wildcard, wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        
        if dlg.ShowModal() == wx.ID_OK:
            filepath = dlg.GetPath()
            
            gridinfo, bodyinfo = readfile(filepath)
            
            self.xygrid = Grid(*gridinfo)
            
            for i, bi in enumerate(bodyinfo):
                self.grid.AppendRows()
                for j, b in enumerate(bi):
                    self.grid.SetCellValue(i, j, to_str(b))
            
            self.grid.ForceRefresh()
        
        dlg.Destroy()

    def on_save(self, event):
        wildcard = "GPrism (*.gprism)|*.gprism"
        dlg = wx.FileDialog(self, "Salvar arquivo", "", "", wildcard, wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        
        if dlg.ShowModal() == wx.ID_OK:
            filepath = dlg.GetPath()
            
            writefile(filepath, self.xygrid, self.prisms)
        
        dlg.Destroy()

    def on_savefig(self, event):
        wildcard = "PNG (*.png)|*.png"
        dlg = wx.FileDialog(self, "Salvar figura", "", "", wildcard, wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        
        if dlg.ShowModal() == wx.ID_OK:
            figurepath = dlg.GetPath()
            
            self.figure.savefig(figurepath, dpi=300)
        
        dlg.Destroy()
    
    def on_export(self, event):
        wildcard = "XYZ (*.xyz)|*.xyz"
        dlg = wx.FileDialog(self, "Exportar dados", "", "", wildcard, wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            filepath = dlg.GetPath()
            
            X, Y, Z = self._getXYZ()
            
            nx, ny = X.shape
            
            with open(filepath, 'w') as f:
                for i in range(nx):
                    for j in range(ny):
                        print >> f, X[i, j], Y[i, j], Z[i, j]
        
        dlg.Destroy()

    def on_add(self, event):
        self.grid.AppendRows()
        self.grid.ForceRefresh()

    def on_remove(self, event):
        to_remove = self.grid.GetSelectedRows()
        self.grid.ClearSelection()

        for i in to_remove[::-1]:
            self.grid.DeleteRows(i)

        self.grid.ForceRefresh()
    
    def on_redraw(self, event):
        self._redraw(self.figure)
    
    def on_config(self, event):
        dlg = ConfigDialog(self)
        dlg.set_grid(self.xygrid)
        fig_props = self._loadfigprops()
        dlg.set_fig_props(fig_props)
        if dlg.ShowModal() == wx.ID_OK:
            self.xygrid = dlg.get_grid()
            fig_props = dlg.get_fig_props()
            self._savefigprops(fig_props)
        
        dlg.Destroy()
    
    def _loadfigprops(self):
        with open(self._figpropspath, 'r') as f:
            figprops = json.load(f)
        
        return figprops
    
    def _savefigprops(self, figprops):
        with open(self._figpropspath, 'w') as f:
            json.dump(figprops, f, sort_keys=True, indent=2, separators=(',', ': '))
    
    def _getXYZ(self):
        x = self.xygrid.xo + np.arange(self.xygrid.nx)*self.xygrid.dx
        y = self.xygrid.yo + np.arange(self.xygrid.ny)*self.xygrid.dy
        X, Y = np.meshgrid(x, y, indexing='ij')
        Z = batchgfield(X, Y, self.prisms)
        
        return X, Y, Z
    
    def _redraw(self, figure):
        X, Y, G = self._getXYZ()
        
        fig_props = self._loadfigprops()
        
        rect = [float(pos)/100.0 for pos in fig_props['pos']]
        rect[2] -= rect[0]
        rect[3] -= rect[1]

        figure.clear()
        ax = figure.add_axes(rect)
        ax.set_aspect('equal', 'box')
        ax.set_title(fig_props['title'])
        ax.set_xlabel(fig_props['xlabel'])
        ax.set_ylabel(fig_props['ylabel'])

        if fig_props['contour']:
            cs = ax.contour(X, Y, G, fig_props['levels'], colors=fig_props['color'])
            if fig_props['clabel']:
                ax.clabel(cs)

        if fig_props['fill']:
            cmap = getattr(matplotlib.cm, fig_props['cmap'], None)
            extent = X[:, 0][0], X[:, 0][-1], Y[0][0], Y[0][-1]
            im = ax.imshow(G.T, origin='lower', extent=extent, cmap=cmap)
            if fig_props['cbar']:
                print "NOT!"

        self.canvas.draw()


if __name__ == '__main__':
    app = wx.App(False)
    mw = MainWindow(None)
    mw.Show()
    app.MainLoop()

