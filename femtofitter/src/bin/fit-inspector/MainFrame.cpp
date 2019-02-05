


#include "Projectionist.hpp"

#include "MainFrame.hpp"

#include <TGStatusBar.h>
#include <TGFileDialog.h>
#include <TGComboBox.h>
#include <TGTextEntry.h>

#include <TG3DLine.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGSlider.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TH3I.h>

#include <TPython.h>
#include <Python.h>

// #include <TGLayoutHints.h>

#include <TString.h>

#include <array>
#include <tuple>
#include <vector>
#include <string>
#include <list>


struct DataQuery {
  std::string cfg;
  std::string cent;
};

struct DataTree {
  // std::map<std::string, // cfg
  //   std::map<std::string, // cent
  //     std::map<std::string, // kt
  //      std::map<std::string, // pair
  //       std::map<std::string, std::string>>>>> root; // mag

  std::vector<std::pair<std::string, // cfg
    std::vector<std::pair<std::string, // cent
      std::vector<std::pair<std::string, // kt
       std::vector<std::pair<std::string, // pair
        std::vector<std::pair<std::string, std::string //magfield
        >> >> >> >> >> root;

  void insert_pieces(const std::vector<std::string> &l)
  {
    auto cfg = l[0],
         cent = l[1],
         kt = l[2],
         pair = l[3],
         magfield = l[4],
         path = l[5];
    // tree.root[cfg][cent][kt][pair][magfield] = path;

    // auto found = std::find(root.begin(), root.end(), [cfg] (auto key) { return key->first == cfg; });
    std::pair<std::string, std::string> mp {magfield, path};

    std::pair<std::string, std::vector<decltype(mp)>> pairp {pair, {mp}};
    std::pair<std::string, std::vector<decltype(pairp)>> ktp {kt, {pairp}};
    std::pair<std::string, std::vector<decltype(ktp)>> centp {cent, {ktp}};
    std::pair<std::string, std::vector<decltype(centp)>> cfg_p {cfg, {centp}};

    // for (auto cfgit = root.begin(); cfgit != root.end(); ++cfgit) {
    for (auto &cfgit : root) {
      if (cfgit.first == cfg) {
      for (auto &cent_it : cfgit.second) {
        if (cent_it.first == cent) {
        for (auto &kt_it : cent_it.second) {
          if (kt_it.first == kt) {
          for (auto &pair_it : kt_it.second) {
            if (pair_it.first == pair) {
            for (auto &mfield_it : pair_it.second) {
              if (mfield_it.second == magfield) {
                std::cerr << "Unexpected duplicate path\n";
                return;
              }
            }
            pair_it.second.emplace_back(mp);
            return;
            }
          }
          kt_it.second.emplace_back(pairp);
          return;
          }
        }
        cent_it.second.emplace_back(ktp);
        return;
        }
      }
      cfgit.second.emplace_back(centp);
      return;
      }
    }
    root.emplace_back(cfg_p);
    return;
  }
};

struct FrameData {
  TString filename;
  TString sfname;

  TGLabel *lbl_fname;
  TGLabel *lbl_sname;

  TGLabel *lbl_xval,
          *lbl_yval,
          *lbl_zval;

  TGHSlider *slider_x,
            *slider_y,
            *slider_z;

  DataTree available_data;

  TGComboBox *cfg_ddown,
             *cent_ddown,
             *kt_ddown,
             *magfield_ddown,
             *pair_ddown;

  void SetFitResultFilename(const TString &fname)
    {
      lbl_fname->SetText(fname.Data());
    }

  void SetSourceFilename(const TString &fname)
    {
      lbl_sname->SetText(fname.Data());
    }

  std::string get_selected_path(const TString &fff)
    {
      const auto
        cfg_id = cfg_ddown->GetSelected(),
        cent_id = cent_ddown->GetSelected(),
        kt_id = kt_ddown->GetSelected(),
        pair_id = pair_ddown->GetSelected(),
        mfld_id = magfield_ddown->GetSelected();

      // const auto
      //   cfgpair = root[cfg_id],
      //   centpair = cfgpair.second[cent_id];

      return root[cfg_id].second[cent_id].second[kt_id].second[pair_id].second[mfld_id].second;
      // return "centpair"
    }

  void build_slider_frame(MyMainFrame *self, TGCompositeFrame *parent)
    {
      TGHorizontalFrame *slidergroup = new TGHorizontalFrame(parent, 300, 50, kHorizontalFrame);
        slidergroup->SetName("slidergroup");

      auto add_labeled_slider = [&self](TGCompositeFrame *p, TString name)
        {
          TGVerticalFrame *frame = new TGVerticalFrame(p,100,50,kVerticalFrame);
           TGHorizontalFrame *lblframe = new TGHorizontalFrame(frame,318,19,kHorizontalFrame);
          //  lblframe->SetName("lblframe");
           TGLabel *axislbl = new TGLabel(lblframe, name);
             axislbl->SetName(Form("lbl%s", name.Data()));
             axislbl->SetTextJustify(36);
             axislbl->SetMargins(0,0,0,0);
             axislbl->SetWidth(100);
             axislbl->SetWrapLength(-1);
             lblframe->AddFrame(axislbl, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,0,2,0));
          //  TGLabel *lblbin = new TGLabel(lblframe, "47");
          //    lblbin->SetTextJustify(36);
          //    lblbin->SetMargins(0,0,0,0);
          //    lblbin->SetWrapLength(-1);
          //    lblframe->AddFrame(lblbin, new TGLayoutHints(kLHintsRight | kLHintsTop,2,2,2,0));
           TGLabel *lblval = new TGLabel(lblframe, "51");
             lblval->SetTextJustify(36);
             lblval->SetMargins(40,50,0,0);
             lblval->SetWrapLength(-1);
             lblframe->AddFrame(lblval, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,0));
           frame->AddFrame(lblframe, new TGLayoutHints(kLHintsCenterX | kLHintsTop | kLHintsExpandX,2,2,2,0));

           TGHSlider *slider = new TGHSlider(frame, 0, kSlider1 | kScaleBoth, -1, kHorizontalFrame);
           slider->SetName(Form("slider_%s", name.Data()));
           slider->SetRange(0, 100);
           slider->SetPosition(51);
           frame->AddFrame(slider, new TGLayoutHints(kLHintsCenterX | kLHintsBottom | kLHintsExpandY| kLHintsExpandX,1,1,0,2));
          // canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","MyMainFrame",this,
          //              "EventInfo(Int_t,Int_t,Int_t,TObject*)");
            slider->Connect("PositionChanged(int)", "MyMainFrame", self, "OnSliderUpdate()");

           p->AddFrame(frame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX));
           return std::make_tuple(slider, lblval);
        };

      std::tie(slider_x, lbl_xval) = add_labeled_slider(slidergroup, "qOut");
      { auto vline = new TGVertical3DLine(slidergroup); slidergroup->AddFrame(vline, new TGLayoutHints(kLHintsExpandY)); }
      std::tie(slider_y, lbl_yval) = add_labeled_slider(slidergroup, "qSide");
      { auto vline = new TGVertical3DLine(slidergroup); slidergroup->AddFrame(vline, new TGLayoutHints(kLHintsExpandY)); }
      std::tie(slider_z, lbl_zval) = add_labeled_slider(slidergroup, "qLong");

      parent->AddFrame(slidergroup, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsTop,2,2,2,2));
      slidergroup->MoveResize(0,0,400,50);
    }

  void build_data_select_box(MyMainFrame *self, TGCompositeFrame *parent)
    {
      auto add_combo_box = [&](int id, TString label)
        {
          TGLabel *lbl = new TGLabel(parent, label);
            lbl->SetTextJustify(36);
            lbl->SetMargins(0,0,8,0);
            lbl->SetWrapLength(-1);
            parent->AddFrame(lbl, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2,2,0,0));

          TGComboBox *cbox = new TGComboBox(parent, id, kHorizontalFrame | kSunkenFrame | kOwnBackground);
            // cbox->SetName("cent_box");
            // fComboBox687->AddEntry("00_05", 0);
            // fComboBox687->AddEntry("05_10", 1);
            cbox->Resize(247,25);
            parent->AddFrame(cbox, new TGLayoutHints(kLHintsNormal, 10));
            cbox->Connect("Selected(int, int)", "MyMainFrame", self, "OnDropdownSelection(int, int)");

          return cbox;
        };

      cfg_ddown = add_combo_box(1, "Config Hash");
      cent_ddown = add_combo_box(2, "Centrality");
      kt_ddown = add_combo_box(3, "kT");
      pair_ddown = add_combo_box(4, "Pair");
      magfield_ddown = add_combo_box(5, "Magnetic Field");
    }

  void build_fitparam_box(MyMainFrame *self, TGCompositeFrame *parent)
    {
      TGLabel *fLabel639 = new TGLabel(parent, "Ro");
        fLabel639->SetTextJustify(36);
        fLabel639->SetMargins(0,0,0,0);
        fLabel639->SetWrapLength(-1);
        parent->AddFrame(fLabel639, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2,2,2,2));
      parent->SetLayoutManager(new TGVerticalLayout(parent));
      // parent->Resize(500, 96);
    }

  void UpdateLabels()
    {
      lbl_xval->SetText(Form("%d  |  %.02f ", slider_x->GetPosition(), float(slider_x->GetPosition()) / slider_x->GetMaxPosition()));
      lbl_yval->SetText(Form("%d  |  %.02f ", slider_y->GetPosition(), float(slider_y->GetPosition()) / slider_y->GetMaxPosition()));
      lbl_zval->SetText(Form("%d  |  %.02f ", slider_z->GetPosition(), float(slider_z->GetPosition()) / slider_z->GetMaxPosition()));
    }

  void update_choices(std::vector<std::string> cfg, std::vector<std::string> cent, std::vector<std::string> kt)
    {
      cent_ddown->RemoveAll();
      for (int i=0; i < cent.size(); ++i) {
        cent_ddown->AddEntry(cent[i].c_str(), i);
      }

      kt_ddown->RemoveAll();
      for (int i=0; i < kt.size(); ++i) {
        kt_ddown->AddEntry(kt[i].c_str(), i);
      }
    }

  void reset_datatree(DataTree tree)
  {
    available_data = tree;

    cfg_ddown->RemoveAll();
    cent_ddown->RemoveAll();
    kt_ddown->RemoveAll();
    pair_ddown->RemoveAll();
    magfield_ddown->RemoveAll();

    int idx_cfg = 0;
    for (auto &cfg_pair : available_data.root) {
      cfg_ddown->AddEntry(cfg_pair.first.c_str(), idx_cfg++);

      if (idx_cfg == 1) {
        cfg_ddown->Select(0, false);

        int idx_cent = 0;
        for (auto &cent_pair : cfg_pair.second) {
          cent_ddown->AddEntry(cent_pair.first.c_str(), idx_cent++);
          if (idx_cent == 1) {
            cent_ddown->Select(0, false);

        int idx_kt = 0;
        for (auto &kt_pair : cent_pair.second) {
          kt_ddown->AddEntry(kt_pair.first.c_str(), idx_kt++);
          if (idx_kt == 1) {
            kt_ddown->Select(0, false);

        int idx_pair = 0;
        for (auto &pair_pair : kt_pair.second) {
          pair_ddown->AddEntry(pair_pair.first.c_str(), idx_pair++);
          if (idx_pair == 1) {
            pair_ddown->Select(0, false);

        int idx_magfield = 0;
        for (auto &magfield_pair : pair_pair.second) {
          magfield_ddown->AddEntry(magfield_pair.first.c_str(), idx_magfield++);
          if (idx_magfield == 1) {
            magfield_ddown->Select(0, true);
          } } } } } } }
        }
      }
    }

  for (auto &ddown : {cfg_ddown,
    cent_ddown,
    kt_ddown,
    pair_ddown,
    magfield_ddown}) {
    std::cout << "/" << ddown->GetSelected();
  }
  std::cout << "\n";

  }
};



ClassImp(MyMainFrame);


MyMainFrame::MyMainFrame(const TGWindow *p)
  : TGMainFrame(p, 600, 400)
  , fCanvas(nullptr)
  , data(new FrameData())
{
  TGVerticalFrame *subframe = new TGVerticalFrame(this, 600, 400);

  TGHorizontalFrame *topbar = new TGHorizontalFrame(subframe, 200, 80);
    TGTextButton *open_button = new TGTextButton(topbar, "&Open");
    open_button->Connect("Clicked()", "MyMainFrame", this, "OnOpen()");
    topbar->AddFrame(open_button, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 50));

    TGVerticalFrame *fname_frame = new TGVerticalFrame(topbar, 400, 80);
      data->lbl_fname = new TGLabel(fname_frame, "<nofile>");
      fname_frame->AddFrame(data->lbl_fname, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));
      data->lbl_fname->SetTextJustify(kTextLeft);

      data->lbl_sname = new TGLabel(fname_frame, "<nofile>");
      data->lbl_sname->SetTextJustify(kTextLeft);
      fname_frame->AddFrame(data->lbl_sname, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));
      topbar->AddFrame(fname_frame, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));

    subframe->AddFrame(topbar, new TGLayoutHints(kLHintsExpandX));

  // { auto line = new TGHorizontal3DLine(topbar); topbar->AddFrame(line, new TGLayoutHints(kLHintsExpandX)); }

  TGHorizontalFrame *mainframe = new TGHorizontalFrame(subframe, 600, 200);
    TGVerticalFrame *leftbar = new TGVerticalFrame(mainframe);
      TGGroupFrame *databox = new TGGroupFrame(leftbar, "Data Select");
      data->build_data_select_box(this, databox);
      leftbar->AddFrame(databox);

      TGGroupFrame *parambox = new TGGroupFrame(leftbar, "Fit Parameters");
      data->build_fitparam_box(this, parambox);
      leftbar->AddFrame(parambox, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX));

      mainframe->AddFrame(leftbar, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY, 2,2,2,2));

    auto *slider_canvas_container = new TGVerticalFrame(mainframe, 400, 200);
      data->build_slider_frame(this, slider_canvas_container);

      fCanvas = new TRootEmbeddedCanvas(0, slider_canvas_container, 500, 400);
      Int_t wid = fCanvas->GetCanvasWindowId();
      auto canvas = new TCanvas("MyCanvas", 10, 10, wid);
      // canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","MyMainFrame",this,
      //              "EventInfo(Int_t,Int_t,Int_t,TObject*)");

      fCanvas->AdoptCanvas(canvas);

      auto canvas_props = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,0,0,1,1);
      slider_canvas_container->AddFrame(fCanvas, canvas_props);

    mainframe->AddFrame(slider_canvas_container, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY | kLHintsExpandX, 2,2,2,2));

  subframe->AddFrame(mainframe, new TGLayoutHints(kLHintsExpandX  | kLHintsExpandY, 0, 0, 0, 0));
  Int_t parts[] = {45, 15, 10, 30};
  auto fStatusBar = new TGStatusBar(subframe, 50, 10, kVerticalFrame);
  fStatusBar->SetParts(parts, 4);
  fStatusBar->Draw3DCorner(kFALSE);
  subframe->AddFrame(fStatusBar, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));

  TGHorizontalFrame *hframe = new TGHorizontalFrame(subframe, 10, 40);

   TGTextButton *draw = new TGTextButton(hframe, "&Draw");
   draw->Connect("Clicked()", "MyMainFrame", this, "OnDraw()");

   hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   subframe->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));


  //  AddFrame(draw, new TGTextButton(hframe, "&Draw "))

  AddFrame(subframe, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY | kLHintsExpandX | kLHintsExpandY));

  SetWindowName("FitInspector");

  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
  Resize(1266, 1051);
}


MyMainFrame::~MyMainFrame()
{
  Cleanup();
  delete fCanvas;
}

void
MyMainFrame::OnDraw()
{
  std::cout << "OnDraw\n";

  TCanvas &c = *fCanvas->GetCanvas();

  c.SetFillColor(42);

  c.Modified();
  c.Update();
}


void
MyMainFrame::OnOpen()
{
  std::cout << "OnOpen\n";

  TGFileInfo fi;
  // fi.fFileTypes = dnd_types;

  // auto filepicker = std::make_unique<TGFileDialog>(
  auto filepicker = new TGFileDialog(
    gClient->GetRoot(),
    this,
    kFDOpen,
    &fi);

  TString selected_filename(fi.fFilename);
  if (selected_filename.IsWhitespace()) {
    return;
  }

  LoadJsonFile(selected_filename);
}


void
MyMainFrame::LoadJsonFile(TString filename)
{
  TPython::Exec("import json");
  TPython::Exec("import pandas as pd");

  TPython::Exec(Form("data = json.load(open('%s'))", filename.Data()));

  std::cout << "Loading filename: " << filename << "\n";

  char *root_file = TPython::Eval("data['filename']");

  TPython::Exec("df = pd.DataFrame(data['df'])");

  data->SetSourceFilename(root_file);
  data->SetFitResultFilename(filename.Data());

  auto load_string_vec = [] (TString cmd)
    {
      std::vector<std::string> result;
      PyObject *pylist = TPython::Eval(cmd);
      auto count = PyList_Size(pylist);
      for (int i=0; i<count; ++i) {
        auto pystr = PyList_GetItem(pylist, i);
        Py_ssize_t len = 0;
        const char *s = PyUnicode_AsUTF8AndSize(pystr, &len);
        result.emplace_back(s, s+len);
      }
      Py_DECREF(pylist);
      return result;
    };

  auto load_string_vec_vec = [&] (TString cmd)
    {
      std::vector<std::vector<std::string>> result;
      PyObject *pylist = TPython::Eval(cmd);
      auto count = PyList_Size(pylist);
      for (int i=0; i<count; ++i) {
        auto pylist_inner = PyList_GetItem(pylist, i);
        auto inner_count = PyList_Size(pylist_inner);

        std::vector<std::string> inresult;
        for (int j=0; j<inner_count; ++j) {
          auto pystr = PyList_GetItem(pylist_inner, j);
          Py_ssize_t len = 0;
          const char *s = PyUnicode_AsUTF8AndSize(pystr, &len);
          inresult.emplace_back(s, s+len);
        }
        result.emplace_back(std::move(inresult));

        Py_DECREF(pylist_inner);
      }
      Py_DECREF(pylist);
      return result;
    };

  std::vector<std::string> cent = load_string_vec("list(df.cent.unique())");
  std::vector<std::string> kt = load_string_vec("list(df.kt.unique())");

  DataTree tree;

  auto vv = load_string_vec_vec("[(list(l) + ['/'.join(l)])"
                                " for l in df[['cfg', 'cent', 'kt', 'pair', 'magfield']].drop_duplicates().values]");

// std::cout << "Loaded: " << vv.size();
  for (auto &l : vv) {

    tree.insert_pieces(l);
  }

  data->reset_datatree(tree);
  std::cout << "Done.\n";
  // data->update_choices(cent, kt);
}

void
MyMainFrame::OnSliderUpdate()
{
  // gSystem->Sleep(100);
  std::cout << "sliders: "
            << data->slider_x->GetPosition() << " "
            << data->slider_y->GetPosition() << " "
            << data->slider_z->GetPosition() << "\n";

  data->UpdateLabels();
}

void
MyMainFrame::OnDropdownSelection(int id, int entry)
{
  std::cout << "DropdownSelection (" << id << ", " << entry << ")\n";

}
